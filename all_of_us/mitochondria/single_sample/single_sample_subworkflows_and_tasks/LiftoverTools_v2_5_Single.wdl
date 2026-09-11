version 1.0

task ChainSwap {
    # a shortcut method (that may be better behaved) which simply swaps reference and target
    input {
        File source_chain
        String input_source_name
        String input_target_name
        Int? preemptible_tries
        String ucsc_docker
    }

    Int disk_size = ceil(size(source_chain, "GB")) * 4
    String chain_output = input_target_name + "_to_" + input_source_name + ".chain"

    command {
        chainSwap ~{source_chain} ~{chain_output}
    }

    output {
        File chain = "~{chain_output}"
    }

    runtime {
        disks: "local-disk " + disk_size + " HDD"
        memory: "1 GB"
        docker: ucsc_docker
        preemptible: select_first([preemptible_tries, 5])
    }
}

task LiftOverVcf {
    # A specialized routine to return the resultant VCF back to GRCh38
    input {
        File original_filtered_vcf
        File new_self_ref_vcf
        File reversed_hom_ref_vcf
        String individual_name

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        File self_fasta
        File self_fasta_index
        File self_dict

        File chain_self_to_ref
        File chain_ref_to_self

        File HailLiftover

        # runtime
        Int? preemptible_tries
        String genomes_cloud_docker
    }

    String this_basename = basename(new_self_ref_vcf, ".vcf")
    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(self_fasta, "GB") + size(self_fasta_index, "GB")
    Int disk_size = ceil(size(original_filtered_vcf, "GB") + size(new_self_ref_vcf, "GB") + size(reversed_hom_ref_vcf, "GB") + ref_size) *2 + 20
    String logging = this_basename + "_fix_liftover.log"

    command <<<
        set -e

        bgzip -c ~{new_self_ref_vcf} > ~{new_self_ref_vcf}.bgz && tabix ~{new_self_ref_vcf}.bgz
        tabix ~{reversed_hom_ref_vcf}
        bcftools isec -p intersected_vcfs -Ov ~{new_self_ref_vcf}.bgz ~{reversed_hom_ref_vcf}

        # there should be no records private to reversed hom ref VCF
        export private_to_rev_hom_ref=$(cat ./intersected_vcfs/0001.vcf | grep ^chrM | wc -l | sed 's/^ *//g')
        if [ $private_to_rev_hom_ref -ne 0 ]; then
          echo "ERROR: There should not be any variants private to the reversed hom ref VCF."
          exit 1;
        fi

        java -jar /usr/gitc/picard.jar LiftoverVcf \
          I=./intersected_vcfs/0000.vcf \
          O=~{this_basename}.selfToRef.pre.vcf \
          R=~{ref_fasta} \
          CHAIN=~{chain_self_to_ref} \
          REJECT=~{this_basename}.selfToRef.rejected.pre.vcf

        java -jar /usr/gitc/picard.jar MergeVcfs \
          I=~{this_basename}.selfToRef.rejected.pre.vcf \
          I=./intersected_vcfs/0002.vcf \
          O=~{this_basename}.selfToRef.rejected.vcf

        sed -e 's/^chr//' ~{original_filtered_vcf} \
          | awk '{OFS="\t"; if (!/^#/){print $1,$2-(length($4)>1 ? 0 : 1),$2-1+length($4),$4"/"$5,"+",length($4),length($5)}}' > filtered.bed
        sed -e 's/^chr//' ~{this_basename}.selfToRef.pre.vcf \
          | awk '{OFS="\t"; if (!/^#/){print $1,$2-(length($4)>1 ? 0 : 1),$2-1+length($4),$4"/"$5,"+",length($4),length($5)}}' > success.bed
        export n_ref_pass_thru=$(bedtools intersect -a filtered.bed -b success.bed | awk '{OFS="\t"; if (($7 == 1)) {print}}' | wc -l | sed 's/^ *//g')
        echo $n_ref_pass_thru > n_ref_pass_thru.txt

        if [ $n_ref_pass_thru -ne 0 ]; then
          echo "ERROR: All variants changed in the self-reference, and all sites within self-reference insertions (excluding the first base), should have failed LiftoverVcf, which is not the case here."
          exit 1;
        fi

        export n_filtered=$(cat ~{original_filtered_vcf} | grep ^chrM | wc -l | sed 's/^ *//g')
        export n_original=$(cat ~{new_self_ref_vcf} | grep ^chrM | wc -l | sed 's/^ *//g')
        export n_pass=$(cat ~{this_basename}.selfToRef.pre.vcf | grep ^chrM | wc -l | sed 's/^ *//g')
        export n_rejected=$(cat ~{this_basename}.selfToRef.rejected.vcf | grep ^chrM | wc -l | sed 's/^ *//g')
        export n_pass_rejected=$((n_pass + n_rejected))
        echo $n_pass > n_pass.txt

        if [ $n_filtered -gt $n_rejected ]; then
          echo "ERROR: the number of sites changed in the self-reference should be less than or equal to the number of sites rejected on initial LiftoverVcf run."
          exit 1;
        fi

        if [ $n_pass_rejected -ne $n_original ]; then
          echo "ERROR: Records appear to have been lost. The sum of records in passing and rejected VCFs should be the same as the number in the original VCF."
          exit 1;
        fi

        echo "$n_filtered sites are changed in self-reference. $n_original sites were found variant after second-round Mutect. Of these, $n_pass passed first-round Liftover and $n_rejected failed and are being piped to Hail pipeline for rescue."

        # now run hail script to fix the rejects
        python3.7 ~{HailLiftover} \
        --vcf-file ~{this_basename}.selfToRef.rejected.vcf \
        --success-vcf-file ~{this_basename}.selfToRef.pre.vcf \
        --self-homoplasmies ~{reversed_hom_ref_vcf} \
        --individual-name ~{individual_name} \
        --self-to-ref-chain ~{chain_self_to_ref} \
        --ref-to-self-chain ~{chain_ref_to_self} \
        --self-fasta ~{self_fasta} \
        --self-fai ~{self_fasta_index} \
        --reference-fasta ~{ref_fasta} \
        --reference-fai ~{ref_fasta_index} \
        --output-prefix ~{this_basename}.round2liftover \
        --export-homoplasmic-deletions-coverage \
        --output-txt-for-wdl \
        --logging ~{logging}

        bgzip -cd ~{this_basename}.round2liftover.rejected.vcf.bgz > ~{this_basename}.round2liftover.rejected.vcf
        bgzip -cd ~{this_basename}.round2liftover.fixed.vcf.bgz > ~{this_basename}.round2liftover.fixed.vcf
        bgzip -cd ~{this_basename}.round2liftover.updated_success.vcf.bgz > ~{this_basename}.round2liftover.updated_success.vcf

        java -jar /usr/gitc/picard.jar MergeVcfs \
          I=~{this_basename}.round2liftover.updated_success.vcf \
          I=~{this_basename}.round2liftover.fixed.vcf \
          O=~{this_basename}.selfToRef.final.vcf

        export n_final_pass=$(cat ~{this_basename}.selfToRef.final.vcf | grep ^chrM | wc -l | sed 's/^ *//g')
        echo $n_final_pass > n_final_pass.txt
    >>>

    runtime {
        disks: "local-disk " + disk_size + " HDD"
        memory: "2500 MB"
        docker: genomes_cloud_docker
        preemptible: select_first([preemptible_tries, 5])
    }

    output {
        File liftover_r1_rejected_vcf = "~{this_basename}.selfToRef.rejected.vcf"
        File liftover_r1_vcf = "~{this_basename}.selfToRef.pre.vcf"
        File liftover_r2_success_r1_vcf = "~{this_basename}.round2liftover.updated_success.vcf"
        File liftover_r2_rejected_vcf = "~{this_basename}.round2liftover.rejected.vcf"
        File liftover_r2_intermediate_vcf = "~{this_basename}.round2liftover.fixed.vcf"
        File liftover_r2_final_vcf = "~{this_basename}.selfToRef.final.vcf"
        File liftover_r2_log = "~{logging}"
        File gap_coverage = "~{this_basename}.round2liftover.deletions_coverage.tsv"
        # stats
        Int n_liftover_changed_selfref_and_passed = read_int('n_ref_pass_thru.txt')
        Int n_liftover_r1_failed = read_int('~{this_basename}.round2liftover.round1_failed_sites.txt')
        Int n_liftover_r1_pass = read_int('n_pass.txt')
        Int n_liftover_r2_failed = read_int('~{this_basename}.round2liftover.round2_failed_sites.txt')
        Int n_liftover_r2_fixed = read_int('~{this_basename}.round2liftover.round2_fixed_sites.txt')
        Int n_liftover_r2_pass = read_int('n_final_pass.txt')
        Int n_liftover_r2_left_shift = read_int('~{this_basename}.round2liftover.left_alignment_of_indels.txt')
        Int n_liftover_r2_injected_from_success = read_int('~{this_basename}.round2liftover.round2_success_injected.txt')
        Int n_liftover_r2_fixed_snv_mismatch = read_int('~{this_basename}.round2liftover.snvs_with_mismatchrefallele.txt')
        Int n_liftover_r2_snv_complex = read_int('~{this_basename}.round2liftover.snvs_complex_liftover.txt')
        Int n_liftover_r2_fixed_ref_insertions = read_int('~{this_basename}.round2liftover.homoplasmic_insertions_rel_to_grch38.txt')
        Int n_liftover_r2_ref_insertion_new_haplo = read_int('~{this_basename}.round2liftover.ref_insertion_new_haplos.txt')
        Int n_liftover_r2_fixed_ref_dele = read_int('~{this_basename}.round2liftover.homoplasmic_deletions_rel_to_grch38.txt')
        Int n_liftover_r2_fixed_indel_first_allele_diff = read_int('~{this_basename}.round2liftover.heteroplasmic_indels_first_site_swapped_to_grch38.txt')
        Int n_liftover_r2_dele_first_allele_internal_diff = read_int('~{this_basename}.round2liftover.deletions_with_first_allele_and_internal_switches.txt')
        Int n_liftover_r2_fixed_dele_internal_diff = read_int('~{this_basename}.round2liftover.heteroplasmic_deletions_with_internal_sites_swapped_to_grch38.txt')
        Int n_liftover_r2_dele_internal_indels = read_int('~{this_basename}.round2liftover.hom_dele_with_indels_reverted.txt')
        Int n_liftover_r2_dele_span_insert = read_int('~{this_basename}.round2liftover.rescued_heteroplasmic_deletions_spanning_homoplasmic_insertions.txt')
        Int n_liftover_r2_dele_span_insert_internal_diff = read_int('~{this_basename}.round2liftover.het_dele_span_insert_repaired_with_complex_rework.txt')
        Int n_liftover_r2_failed_het_dele_span_insertion_boundary = read_int('~{this_basename}.round2liftover.het_deletions_span_insertions.txt')
        Int n_liftover_r2_failed_new_dupes_leftshift = read_int('~{this_basename}.round2liftover.new_dupes_left_shift_failed.txt')
        Int n_liftover_r2_fancy_flip_af = read_int('~{this_basename}.round2liftover.sites_with_fancy_flip.txt')
        Int n_liftover_r2_repaired_success = read_int('~{this_basename}.round2liftover.success_sites_flipped.txt')
        Int n_liftover_r2_fixed_indelstraddlesintervals_flag = read_int('~{this_basename}.round2liftover.fixed_indelStraddlesMultipleIntevals.txt')
    }
}

task LiftOverCoverage {
    # run liftover provided a coverage txt file
    # use reference allele depth from the SELF-REFERENCE data after reversing the force-call
    # to impute the appropriate coverage for a deleted site. This REF allele is the
    # actual original sequence.
    input {
        String sample_name
        File self_ref_table
        File chain
        File homoplasmic_deletions_coverage
        Int? preemptible_tries
        String ucsc_docker
    }

    String this_basename = sample_name + "_" + basename(self_ref_table, '.tsv')
    String output_bed = this_basename + ".liftedOver.bed"
    String rejects = this_basename + ".liftedOver.unmapped"
    String final_tsv = this_basename + ".appended.liftedOver.tsv"
    Int disk_size = ceil(size(self_ref_table, "GB") * 4) + 20

    command <<<
        tail -n +2 ~{self_ref_table} | awk '{OFS="\t"} {print $1"\t"$2-1"\t"$2"\t"$4}' > per_base_coverage.bed
        liftOver per_base_coverage.bed ~{chain} ~{output_bed} ~{rejects}

        R --vanilla <<CODE
        init_chain <- readLines("~{chain}")
        if (length(grep('^chain', init_chain)) != 1) {
          stop('ERROR: only single block chain files are currently supported.')
        }
        loc_init_chain <- grep('^chain', init_chain)
        chain_splits = strsplit(init_chain[loc_init_chain],'\\\\s+')
        if (!all(sapply(chain_splits, function(x)as.numeric(c(x[6], x[11]))) == 0)) {
          stop('ERROR: chain file must start at 0 for both source and target.')
        }
        if (!(1 %in% loc_init_chain)) {
          stop('ERROR: first line of chain file must be header.')
        }
        len_tab <- length(strsplit(init_chain[2],'\t')[[1]])
        len_spc <- length(strsplit(init_chain[2],' ')[[1]])
        if (len_tab > len_spc) {
          sep <- '\t'
        } else {
          sep <- ' '
        }

        chain <- read.csv("~{chain}", skip=1, header=F, col.names=c('bases_in_block', 'ns', 'nt'), sep=sep, stringsAsFactors=F, fill=T)
        chain[nrow(chain),2:3] <- 0
        if (nrow(chain) == 1) {
          chain[['coord_s']] = chain[['bases_in_block']]
          chain[['coord_t']] = chain[['bases_in_block']]
        } else {
          chain[['coord_s']] = cumsum(chain[['bases_in_block']]) + cumsum(c(0,chain[['ns']][2:nrow(chain)-1]))
          chain[['coord_t']] = cumsum(chain[['bases_in_block']]) + cumsum(c(0,chain[['nt']][2:nrow(chain)-1]))
        }

        if (any(is.na(chain))) {
          stop('ERROR: there should be no missing fields in chain.')
        }
        if (sum((chain[['ns']] > 0) & (chain[['nt']] > 0)) != 0) {
          stop('ERROR: Chain files should not have breaks -- e.g., each ungapped block should be followed by a gap in EITHER the source or target, but not both.')
        }

        start_end_to_bed <- function(start, end_inclusive) {
          # push start up by one, because INDELs share the first base
          new_start = start+1
          ends = end_inclusive
          if (length(new_start) != length(ends)) {
            stop('ERROR: starts and ends must have the same length.')
          }
          if (length(new_start) == 0) {
            empty_df <- data.frame('chr'='chrM', 'start'=1, 'end'=1)
            return(empty_df[empty_df[['start']] > 1,])
          }
          bp = unlist(lapply(1:length(start), function(idx) new_start[idx]:ends[idx]))
          df = data.frame('chr' = 'chrM', 'start' = bp-1, 'end' = bp)
          return(df[order(df[['start']]),])
        }

        bed_entries_missing_in_db <- function(query, db) {
          if ((nrow(query) > 0) & (nrow(db) > 0)) {
            return(query[!paste0(query[['start']],':',query[['end']]) %in% paste0(db[['start']],':',db[['end']]),])
          } else {
            return(query)
          }
        }

        n_rows_different_keys_bed <- function(df1, df2) {
          return(nrow(bed_entries_missing_in_db(df1, df2)) + nrow(bed_entries_missing_in_db(df2, df1)))
        }

        insertions <- chain[chain[['ns']] > 0,]
        expected_failed_liftover <- start_end_to_bed(insertions[['coord_s']], insertions[['coord_s']]+insertions[['ns']])
        failure_bed <- read.csv("~{rejects}", header=F, sep='\t', col.names=c('chr','start','end','cov'), comment.char='#')
        if (n_rows_different_keys_bed(expected_failed_liftover, failure_bed) > 0) {
          stop('ERROR: Rejected bed file is not identical to expectation.')
        }

        target_len <- as.numeric(strsplit(init_chain[1], '\\\\s+')[[1]][9])
        deletions <- chain[chain[['nt']] > 0,]
        expected_liftover_gaps <- start_end_to_bed(deletions[['coord_t']], deletions[['coord_t']]+deletions[['nt']])
        success_bed <- read.csv("~{output_bed}", header=F, sep='\t', col.names=c('chr','start','end','cov'), comment.char='#')
        full_bed <- data.frame('chr'='chrM', 'start'=0:(target_len-1), 'end'=1:target_len)
        not_in_success_bed <- bed_entries_missing_in_db(full_bed, success_bed)
        if (n_rows_different_keys_bed(expected_liftover_gaps, not_in_success_bed) > 0) {
          stop('ERROR: Success bed file has different gaps than expected based on chain.')
        }
        if (nrow(expected_liftover_gaps) == 0) {
          print('NOTE: No gaps found in coverage file, so no changes were made.')
          appended_bed <- success_bed
          appended_bed <- appended_bed[order(appended_bed[['start']]),]
        } else {
          fill_in_coverage <- read.csv("~{homoplasmic_deletions_coverage}", sep='\t', stringsAsFactors=F)
          bed_homoplas_dele <- data.frame('chr'='chrM', 'start'=fill_in_coverage[['reference_position']]-1,
                                          'end'=fill_in_coverage[['reference_position']], 'cov'=fill_in_coverage[['ref_allele_depth']])
          if (n_rows_different_keys_bed(expected_liftover_gaps, bed_homoplas_dele) > 0) {
            stop('ERROR: Records in homoplasmic_deletions_coverage are different from expected based on chain.')
          }
          appended_bed <- rbind(success_bed, merge(not_in_success_bed, bed_homoplas_dele, by=c('chr','start','end')))
          appended_bed <- appended_bed[order(appended_bed[['start']]),]
          row.names(appended_bed) <- NULL
          if ((nrow(bed_entries_missing_in_db(full_bed, appended_bed)) > 0) | (!all(appended_bed[['end']] == 1:target_len))) {
            stop('ERROR: the final lifted over bed file should have one row per position on the mtDNA.')
          }
        }

        final_tsv <- appended_bed
        names(final_tsv) <- c('chrom', 'to_rm', 'pos', 'coverage')
        final_tsv[['target']] <- '.'
        write.table(final_tsv[,c('chrom','pos','target','coverage')], "~{final_tsv}", row.names=F, sep='\t', quote=F)
        CODE
    >>>

    output {
        File reference_coverage = "~{final_tsv}"
        File intermediate_bed = "~{output_bed}"
        File rejected = "~{rejects}"
    }

    runtime {
        disks: "local-disk " + disk_size + " HDD"
        memory: "2 GB"
        docker: ucsc_docker
        preemptible: select_first([preemptible_tries, 5])
    }
}

task LiftOverAndJoinCoverage {
    # As this is meant for QC, currently this does not allow spiked in
    # coverage at sites with gaps
    input {
        File ref_table
        File self_table
        File self_table_shifted
        File? chain
        Int? preemptible_tries
        Int? mem
        Boolean? use_fast_join = true
        String ucsc_docker
    }

    String this_basename = basename(self_table, '.tsv')
    String output_bed = this_basename + ".nuc.liftedOver.bed"
    String output_bed_shifted = this_basename + ".shifted.nuc.liftedOver.bed"
    String rejects = this_basename + ".nuc.liftedOver.unmapped"
    String rejects_shifted = this_basename + ".shifted.nuc.liftedOver.unmapped"
    String final_tsv = this_basename + ".nuc.appended.liftedOver.tsv"
    Int disk_size = ceil(size(self_table, "GB") * 4) + 20

    # for an attempted speed up, we removed the following from below:
    # joint_bed <- merge(ref_pre_bed$tab, self_lifted_bed$tab, by=c('chr', 'start', 'end'), all=TRUE)
    # joint_bed <- merge(joint_bed, self_lifted_shifted_bed$tab, by=c('chr', 'start', 'end'), all=TRUE)

    command <<<
        tail -n +2 ~{self_table} | awk '{OFS="\t"} {print $1"\t"$2-1"\t"$2"\t"$4"\t"$3}' > per_base_coverage.bed
        tail -n +2 ~{self_table_shifted} | awk '{OFS="\t"} {print $1"\t"$2-1"\t"$2"\t"$4"\t"$3}' > per_base_coverage_shifted.bed
        tail -n +2 ~{ref_table} | awk '{OFS="\t"} {print $1"\t"$2-1"\t"$2"\t"$4"\t"$3}' > per_base_coverage_pre.bed
        liftOver per_base_coverage.bed ~{chain} ~{output_bed} ~{rejects}
        liftOver per_base_coverage_shifted.bed ~{chain} ~{output_bed_shifted} ~{rejects_shifted}

        R --vanilla <<CODE
          use_fast_join <- "~{use_fast_join}" == "true"
          ref_pre_bed <- read.csv("per_base_coverage_pre.bed", header=F, sep='\t', col.names=c('chr','start','end','cov_ref','target'), comment.char='#', stringsAsFactors=F)
          self_lifted_bed <- read.csv("~{output_bed}", header=F, sep='\t', col.names=c('chr','start','end','cov_self_ref','target'), comment.char='#', stringsAsFactors=F)
          self_lifted_shifted_bed <- read.csv("~{output_bed_shifted}", header=F, sep='\t', col.names=c('chr','start','end','cov_self_ref_shifted','target'), comment.char='#', stringsAsFactors=F)
          original_cov <- read.csv('~{ref_table}', sep='\t', stringsAsFactors=F)

          edit_self_bed <- function(self_lifted_bed) {
            chr_to_convert <- unique(self_lifted_bed[['chr']])
            mapped_chrs <- lapply(chr_to_convert, function(x) unique(original_cov[['chrom']][original_cov[['target']] == x]))
            if(!all(sapply(mapped_chrs,length) == 1)) {
              stop('ERROR: some intervals did not map to chromosomes in the original reference coverage files (or mapped to multiple).')
            }
            chr_mapper <- sapply(mapped_chrs,function(x)x[[1]])
            names(chr_mapper) <- chr_to_convert
            self_lifted_bed[['chr']] <- chr_mapper[self_lifted_bed[['chr']]]
            split_by_targ <- split(self_lifted_bed, self_lifted_bed[['target']])
            split_by_targ_updated <- lapply(names(split_by_targ), function(x) {
              this_tab <- split_by_targ[[x]]
              filtered_ref_start <- ref_pre_bed[ref_pre_bed[['target']] == x,][['start']]
              this_tab[['start']] <- this_tab[['start']] + filtered_ref_start[1]
              this_tab[['end']] <- this_tab[['end']] + filtered_ref_start[1]
              return(this_tab)
            })
            self_lifted_bed <- do.call("rbind", split_by_targ_updated)
            return(self_lifted_bed)
          }

          ensure_unique <- function(tab, keep_cols=F) {
            tf_redundant <- all(tab[['end']] - tab[['start']]) == 1
            if (tf_redundant) {
              tab[['joiner']] <- paste0(tab[['chr']], ':', tab[['start']], ':', tab[['target']])
              vals <- c('chr', 'start', 'target')
            } else {
              stop('ERROR: end is not the same as start+1 -- this is not handled.')
              tab[['joiner']] <- paste0(tab[['chr']], ':', tab[['start']], ':', tab[['end']], ':', self_lifted_bed[['target']])
              vals <- c('chr', 'start', 'end', 'target')
            }
            if (!keep_cols) {
              tab <- tab[,-1*which(names(tab) %in% vals)]
            }
            if (any(duplicated(tab[['joiner']]))) {
              stop('Joining rows not unique. This means you have duplicate loci in the per-base coverage results.')
            } else {
              return(list('vals' = vals, 'tab' = tab))
            }
          }

          fast_join <- function(tab1, tab2, by, target) {
            all_keys <- unique(c(tab1[[by]], tab2[[by]]))
            new_keys <- all_keys[!all_keys %in% tab1[[by]]]
            if (length(new_keys) > 0) {
              new_df <- data.frame(ff = new_keys)
              names(new_df) <- by
              new_df[setdiff(names(tab1), names(new_df))] <- NA
              tab1 <- rbind(tab1, new_df)
            }
            tab1[[target]] <- tab2[[target]][match(tab1[[by]], tab2[[by]])]
            return(tab1)
          }

          self_lifted_bed <- edit_self_bed(self_lifted_bed)
          self_lifted_shifted_bed <- edit_self_bed(self_lifted_shifted_bed)

          if (use_fast_join) {
            ref_pre_bed <- ensure_unique(ref_pre_bed, T)
            self_lifted_bed <- ensure_unique(self_lifted_bed)
            self_lifted_shifted_bed <- ensure_unique(self_lifted_shifted_bed)
            joint_bed <- fast_join(ref_pre_bed[['tab']], self_lifted_bed[['tab']], by='joiner', target='cov_self_ref')
            joint_bed <- fast_join(joint_bed, self_lifted_shifted_bed[['tab']], by='joiner', target='cov_self_ref_shifted')
            names(joint_bed) <- c('chrom', 'to_rm', 'pos', 'coverage_original', 'target', 'to_rm2', 'coverage_remapped_self', 'coverage_remapped_self_shifted')
          } else {
            joint_bed <- merge(ref_pre_bed, self_lifted_bed, by=c('chr', 'start', 'end', 'target'), all=TRUE)
            joint_bed <- merge(joint_bed, self_lifted_shifted_bed, by=c('chr', 'start', 'end', 'target'), all=TRUE)
            names(joint_bed) <- c('chrom', 'to_rm', 'pos', 'target', 'coverage_original', 'coverage_remapped_self', 'coverage_remapped_self_shifted')
          }

          write.table(joint_bed[,c('chrom','pos','target','coverage_original','coverage_remapped_self','coverage_remapped_self_shifted')], "cov.tsv", row.names=F, sep='\t', quote=F)
        CODE

        bgzip -c cov.tsv > ~{final_tsv}.bgz
    >>>

    output {
        File reference_coverage = "~{final_tsv}.bgz"
        File intermediate_bed = "~{output_bed}"
        File rejected = "~{rejects}"
        File intermediate_bed_shifted = "~{output_bed_shifted}"
        File rejected_shifted = "~{rejects_shifted}"
    }

    runtime {
        disks: "local-disk " + disk_size + " HDD"
        memory: select_first([mem, 500]) + " MB"
        docker: ucsc_docker
        preemptible: select_first([preemptible_tries, 5])
    }
}