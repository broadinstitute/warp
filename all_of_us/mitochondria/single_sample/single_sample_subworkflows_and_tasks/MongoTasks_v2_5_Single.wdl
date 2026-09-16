version 1.0
# This is a stupid file containing giant tasks for use with pipeline platforms
# that cannot handle a lot of parallel complexity. This is obviously hard to maintain.
# These tasks are serial combinations of more modular tasks that were located in each individual WDL.

task MongoSubsetBamToChrMAndRevert {
    input {
        File input_bam
        File? input_bai
        String sample_name

        File? mt_interval_list
        File? nuc_interval_list
        String? contig_name
        String? requester_pays_project
        File? ref_fasta
        File? ref_fasta_index
        File? ref_dict

        Boolean skip_restore_hardclips
        String? read_name_regex
        Int? read_length
        Int? coverage_cap

        File? gatk_override
        String? gatk_docker_override
        String gatk_version
        String? printreads_extra_args

        # runtime
        Boolean force_manual_download # will download using gsutil cp
        Int? mem
        Int? n_cpu
        Int? preemptible_tries
    }
    Float ref_size = if defined(ref_fasta) then size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB") else 0
    String appended_bam = input_bam + '.crai'
    Int disk_size = ceil(ref_size) + ceil(size(input_bam,'GB')) + 20
    Int read_length_for_optimization = select_first([read_length, 151])
    Int machine_mem = select_first([mem, 4])
    Int command_mem = (machine_mem * 1000) - 500
    String skip_hardclip_str = if skip_restore_hardclips then "--RESTORE_HARDCLIPS false" else ""
    String requester_pays_prefix = (if defined(requester_pays_project) then "-u " else "") + select_first([requester_pays_project, ""])

    String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

    meta {
        description: "Subsets a whole genome bam to just Mitochondria reads"
    }
    parameter_meta {
        ref_fasta: "Reference is only required for cram input. If it is provided ref_fasta_index and ref_dict are also required."
        input_bam: {
                       localization_optional: true
                   }
        input_bai: {
                       localization_optional: true
                   }
    }
    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        mkdir out

        this_bam="~{input_bam}"
        this_bai="~{select_first([input_bai, appended_bam])}"
        this_sample=out/"~{sample_name}"

        ~{if force_manual_download then "gsutil " + requester_pays_prefix + " cp ~{d}{this_bam} bamfile.cram" else ""}
        ~{if force_manual_download then "gsutil " + requester_pays_prefix + " cp ~{d}{this_bai} bamfile.cram.crai" else ""}
        ~{if force_manual_download then "this_bam=bamfile.cram" else ""}
        ~{if force_manual_download then "this_bai=bamfile.cram.crai" else ""}

        gatk --java-options "-Xmx~{command_mem}m" PrintReads \
          ~{"-R " + ref_fasta} \
          ~{"-L " + mt_interval_list} \
          ~{"-L " + nuc_interval_list} \
          ~{"-L " + contig_name} \
          --read-filter MateOnSameContigOrNoMappedMateReadFilter \
          --read-filter MateUnmappedAndUnmappedReadFilter \
          ~{"--gcs-project-for-requester-pays " + requester_pays_project} \
          ~{if force_manual_download then '-I bamfile.cram --read-index bamfile.cram.crai' else "-I ~{d}{this_bam} --read-index ~{d}{this_bai}"} \
          -O "~{d}{this_sample}.bam"  ~{printreads_extra_args}

        echo "Now removing mapping..."
        set +e
        gatk --java-options "-Xmx~{command_mem}m" ValidateSamFile \
          -INPUT "~{d}{this_sample}.bam" \
          -O output.txt \
          -M VERBOSE \
          -IGNORE_WARNINGS true \
          -MAX_OUTPUT 9999999
        cat output.txt | \
          grep 'ERROR.*Mate not found for paired read' | \
          sed -e 's/ERROR::MATE_NOT_FOUND:Read name //g' | \
          sed -e 's/, Mate not found for paired read//g' > read_list.txt
        cat read_list.txt | wc -l | sed 's/^ *//g' > "~{d}{this_sample}.ct_failed.txt"
        if [[ $(tr -d "\r\n" < read_list.txt|wc -c) -eq 0 ]]; then
          cp "~{d}{this_sample}.bam" rescued.bam
        else
          gatk --java-options "-Xmx~{command_mem}m" FilterSamReads \
            -I "~{d}{this_sample}.bam" \
            -O rescued.bam \
            -READ_LIST_FILE read_list.txt \
            -FILTER excludeReadList
        fi
        gatk --java-options "-Xmx~{command_mem}m" RevertSam \
          -INPUT rescued.bam \
          -OUTPUT_BY_READGROUP false \
          -OUTPUT "~{d}{this_sample}.unmap.bam" \
          -VALIDATION_STRINGENCY LENIENT \
          -ATTRIBUTE_TO_CLEAR FT \
          -ATTRIBUTE_TO_CLEAR CO \
          -SORT_ORDER queryname \
          -RESTORE_ORIGINAL_QUALITIES false ~{skip_hardclip_str}

        set -e
        echo "Now getting WGS metrics on the subsetted bam..."
        gatk --java-options "-Xmx~{command_mem}m" CollectWgsMetrics \
          INPUT="~{d}{this_sample}.bam" \
          ~{"INTERVALS=" + mt_interval_list} \
          VALIDATION_STRINGENCY=SILENT \
          REFERENCE_SEQUENCE=~{ref_fasta} \
          OUTPUT="~{d}{this_sample}.wgs_metrics.txt" \
          USE_FAST_ALGORITHM=true \
          READ_LENGTH=~{read_length_for_optimization} \
          ~{"COVERAGE_CAP=" + coverage_cap} \
          INCLUDE_BQ_HISTOGRAM=true \
          THEORETICAL_SENSITIVITY_OUTPUT="~{d}{this_sample}.theoretical_sensitivity.txt"

        R --vanilla <<CODE
          df = read.table("~{d}{this_sample}.wgs_metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
          write.table(floor(df[,"MEAN_COVERAGE"]), "~{d}{this_sample}.mean_coverage.txt", quote=F, col.names=F, row.names=F)
          write.table(df[,"MEDIAN_COVERAGE"], "~{d}{this_sample}.median_coverage.txt", quote=F, col.names=F, row.names=F)
        CODE

        echo "Now preprocessing subsetted bam..."
        gatk --java-options "-Xmx~{command_mem}m" MarkDuplicates \
          INPUT="~{d}{this_sample}.bam" \
          OUTPUT=md.bam \
          METRICS_FILE="~{d}{this_sample}.duplicate.metrics" \
          VALIDATION_STRINGENCY=SILENT \
          ~{"READ_NAME_REGEX=" + read_name_regex} \
          OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
          ASSUME_SORT_ORDER="queryname" \
          CLEAR_DT="false" \
          ADD_PG_TAG_TO_READS=false

        gatk --java-options "-Xmx~{command_mem}m" SortSam \
          INPUT=md.bam \
          OUTPUT="~{d}{this_sample}.proc.bam" \
          SORT_ORDER="coordinate" \
          CREATE_INDEX=true \
          MAX_RECORDS_IN_RAM=300000
    >>>
    runtime {
        memory: machine_mem + " GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
        preemptible: select_first([preemptible_tries, 5])
        cpu: select_first([n_cpu,1])
    }
    output {
        File output_bam = "out/~{sample_name}.proc.bam"
        File output_bai = "out/~{sample_name}.proc.bai"
        File unmapped_bam = "out/~{sample_name}.unmap.bam"
        File duplicate_metrics = "out/~{sample_name}.duplicate.metrics"
        Int reads_dropped = read_int("out/~{sample_name}.ct_failed.txt")
        Int mean_coverage = read_int("out/~{sample_name}.mean_coverage.txt")
    }
}

task MongoSubsetBamToChrMAndRevertFUSE {
    input {
        File input_bam
        File? input_bai
        String sample_name

        File? mt_interval_list
        File? nuc_interval_list
        String? contig_name
        String? requester_pays_project
        File? ref_fasta
        File? ref_fasta_index
        File? ref_dict

        Boolean skip_restore_hardclips
        String? read_name_regex
        Int? read_length
        Int? coverage_cap

        File? gatk_override
        String? gatk_docker_override
        String gatk_version

        # runtime
        Boolean force_manual_download # will download using gsutil cp
        Int? mem
        Int? n_cpu
        Int? preemptible_tries
    }
    Float ref_size = if defined(ref_fasta) then size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB") else 0
    Int disk_size = ceil(ref_size) + ceil(size(input_bam,'GB')) + 20
    Int read_length_for_optimization = select_first([read_length, 151])
    Int machine_mem = select_first([mem, 4])
    Int command_mem = (machine_mem * 1000) - 500
    String skip_hardclip_str = if skip_restore_hardclips then "--RESTORE_HARDCLIPS false" else ""
    String requester_pays_prefix = (if defined(requester_pays_project) then "-u " else "") + select_first([requester_pays_project, ""])

    String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

    meta {
        description: "Subsets a whole genome bam to just Mitochondria reads"
    }
    parameter_meta {
        ref_fasta: "Reference is only required for cram input. If it is provided ref_fasta_index and ref_dict are also required."
        input_bam: {
                       localization_optional: true
                   }
        input_bai: {
                       localization_optional: true
                   }
    }
    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        ls -l /
        ls -l /mnt

        mkdir out

        this_bam="~{input_bam}"
        this_bai="~{select_first([input_bai, input_bam + '.crai'])}"
        this_sample=out/"~{sample_name}"

        ~{if force_manual_download then "gsutil " + requester_pays_prefix + " cp ~{d}{this_bam} bamfile.cram" else ""}
        ~{if force_manual_download then "gsutil " + requester_pays_prefix + " cp ~{d}{this_bai} bamfile.cram.crai" else ""}
        ~{if force_manual_download then "this_bam=bamfile.cram" else ""}
        ~{if force_manual_download then "this_bai=bamfile.cram.crai" else ""}

        gatk --java-options "-Xmx~{command_mem}m" PrintReads \
          ~{"-R " + ref_fasta} \
          ~{"-L " + mt_interval_list} \
          ~{"-L " + nuc_interval_list} \
          ~{"-L " + contig_name} \
          --read-filter MateOnSameContigOrNoMappedMateReadFilter \
          --read-filter MateUnmappedAndUnmappedReadFilter \
          ~{"--gcs-project-for-requester-pays " + requester_pays_project} \
          ~{if force_manual_download then '-I bamfile.cram --read-index bamfile.cram.crai' else "-I ~{d}{this_bam} --read-index ~{d}{this_bai}"} \
          -O "~{d}{this_sample}.bam"

        echo "Now removing mapping..."
        set +e
        gatk --java-options "-Xmx~{command_mem}m" ValidateSamFile \
          -INPUT "~{d}{this_sample}.bam" \
          -O output.txt \
          -M VERBOSE \
          -IGNORE_WARNINGS true \
          -MAX_OUTPUT 9999999
        cat output.txt | \
          grep 'ERROR.*Mate not found for paired read' | \
          sed -e 's/ERROR::MATE_NOT_FOUND:Read name //g' | \
          sed -e 's/, Mate not found for paired read//g' > read_list.txt
        cat read_list.txt | wc -l | sed 's/^ *//g' > "~{d}{this_sample}.ct_failed.txt"
        if [[ $(tr -d "\r\n" < read_list.txt|wc -c) -eq 0 ]]; then
          cp "~{d}{this_sample}.bam" rescued.bam
        else
          gatk --java-options "-Xmx~{command_mem}m" FilterSamReads \
            -I "~{d}{this_sample}.bam" \
            -O rescued.bam \
            -READ_LIST_FILE read_list.txt \
            -FILTER excludeReadList
        fi
        gatk --java-options "-Xmx~{command_mem}m" RevertSam \
          -INPUT rescued.bam \
          -OUTPUT_BY_READGROUP false \
          -OUTPUT "~{d}{this_sample}.unmap.bam" \
          -VALIDATION_STRINGENCY LENIENT \
          -ATTRIBUTE_TO_CLEAR FT \
          -ATTRIBUTE_TO_CLEAR CO \
          -SORT_ORDER queryname \
          -RESTORE_ORIGINAL_QUALITIES false ~{skip_hardclip_str}

        set -e
        echo "Now getting WGS metrics on the subsetted bam..."
        gatk --java-options "-Xmx~{command_mem}m" CollectWgsMetrics \
          INPUT="~{d}{this_sample}.bam" \
          ~{"INTERVALS=" + mt_interval_list} \
          VALIDATION_STRINGENCY=SILENT \
          REFERENCE_SEQUENCE=~{ref_fasta} \
          OUTPUT="~{d}{this_sample}.wgs_metrics.txt" \
          USE_FAST_ALGORITHM=true \
          READ_LENGTH=~{read_length_for_optimization} \
          ~{"COVERAGE_CAP=" + coverage_cap} \
          INCLUDE_BQ_HISTOGRAM=true \
          THEORETICAL_SENSITIVITY_OUTPUT="~{d}{this_sample}.theoretical_sensitivity.txt"

        R --vanilla <<CODE
          df = read.table("~{d}{this_sample}.wgs_metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
          write.table(floor(df[,"MEAN_COVERAGE"]), "~{d}{this_sample}.mean_coverage.txt", quote=F, col.names=F, row.names=F)
          write.table(df[,"MEDIAN_COVERAGE"], "~{d}{this_sample}.median_coverage.txt", quote=F, col.names=F, row.names=F)
        CODE

        echo "Now preprocessing subsetted bam..."
        gatk --java-options "-Xmx~{command_mem}m" MarkDuplicates \
          INPUT="~{d}{this_sample}.bam" \
          OUTPUT=md.bam \
          METRICS_FILE="~{d}{this_sample}.duplicate.metrics" \
          VALIDATION_STRINGENCY=SILENT \
          ~{"READ_NAME_REGEX=" + read_name_regex} \
          OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
          ASSUME_SORT_ORDER="queryname" \
          CLEAR_DT="false" \
          ADD_PG_TAG_TO_READS=false

        gatk --java-options "-Xmx~{command_mem}m" SortSam \
          INPUT=md.bam \
          OUTPUT="~{d}{this_sample}.proc.bam" \
          SORT_ORDER="coordinate" \
          CREATE_INDEX=true \
          MAX_RECORDS_IN_RAM=300000
    >>>
    runtime {
        memory: machine_mem + " GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
        preemptible: select_first([preemptible_tries, 5])
        cpu: select_first([n_cpu,1])
    }
    output {
        File output_bam = "out/~{sample_name}.proc.bam"
        File output_bai = "out/~{sample_name}.proc.bai"
        File unmapped_bam = "out/~{sample_name}.unmap.bam"
        File duplicate_metrics = "out/~{sample_name}.duplicate.metrics"
        Int reads_dropped = read_int("out/~{sample_name}.ct_failed.txt")
        Int mean_coverage = read_int("out/~{sample_name}.mean_coverage.txt")
    }
}

task MongoProduceSelfReference {
    # updated 220616 to output a warning for overlapping mtDNA variants and just output the count, along
    # with a filtering pipeline to remove duplicates and prevent downstream errors
    input {
        String sample_name
        File input_nuc_vcf
        File input_mt_vcf

        String suffix
        File ref_fasta
        File ref_fasta_index
        File nuc_interval_list
        File mt_ref_fasta
        File mt_ref_fasta_index
        File mt_interval_list
        File non_control_region_interval_list

        File fa_renaming_script
        File variant_bounds_script
        File check_hom_overlap_script
        Int? preemptible_tries

        Int n_shift

        String genomes_cloud_docker
        String intertext
    }

    String mt_intervals_basename = basename(mt_interval_list, ".interval_list")
    String noncontrol_basename = basename(non_control_region_interval_list, ".interval_list")
    String selfref_intervals_suffix = ".SelfRefLiftover.interval_list"
    String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

    command <<<
        set -e

        mkdir out

        java -jar /usr/gitc/picard.jar IntervalListTools SORT=true I=~{mt_interval_list} O=internal_mt.interval_list
        java -jar /usr/gitc/picard.jar IntervalListTools SORT=true I=~{nuc_interval_list} O=internal_nuc.interval_list

        java -jar /usr/gitc/picard.jar ExtractSequences \
          INTERVAL_LIST=internal_mt.interval_list \
          R=~{mt_ref_fasta} \
          O=mt_fasta.fasta
        Rscript --vanilla ~{fa_renaming_script} mt_fasta.fasta internal_mt.interval_list FALSE mt_fasta_renamed.fasta TRUE
        samtools faidx mt_fasta_renamed.fasta

        this_nuc_vcf="~{input_nuc_vcf}"
        this_mt_vcf="~{input_mt_vcf}"
        this_samp_name="~{sample_name}"
        this_sample=out/"~{sample_name}"

        this_mt_chain=out/reference_to_"~{d}{this_samp_name}.chain"
        this_nuc_chain=out/reference_to_"~{d}{this_samp_name}NucOnly.chain"
        this_basename="~{d}{this_sample}~{suffix}"
        this_basename_both="~{d}{this_sample}andNuc~{suffix}"
        this_basename_nuc="~{d}{this_sample}NucOnly~{suffix}"
        this_shifted_basename="~{d}{this_basename}.shifted_by_~{n_shift}_bases"
        this_shifted_basename_append="~{d}{this_shifted_basename}.cat"

        this_filt_vcf="~{d}{this_sample}.mtref.homoplasmies.vcf"
        this_mt_fasta="~{d}{this_basename}.fasta"
        this_nuc_mt_fasta="~{d}{this_basename_both}.fasta"
        this_nuc_only_fasta="~{d}{this_basename_nuc}.fasta"
        this_fasta_shifted="~{d}{this_shifted_basename}.fasta"
        this_fasta_cat_shifted="~{d}{this_shifted_basename_append}.fasta"
        this_chain_shifted="~{d}{this_shifted_basename}.shift_back_~{n_shift}_bases.chain"
        this_chain_fwd_shifted="~{d}{this_shifted_basename}.shift_fwd_~{n_shift}_bases.chain"
        this_mt_intervals=out/"~{mt_intervals_basename}.~{d}{this_samp_name}~{selfref_intervals_suffix}"
        this_nonctrl_interval=out/"~{noncontrol_basename}.~{d}{this_samp_name}~{selfref_intervals_suffix}"
        this_ctrl_interval=out/control_region_shifted.chrM."~{d}{this_samp_name}~{selfref_intervals_suffix}"

        this_vcf_bgz="~{d}{this_basename}.reversed.selfRef.homoplasmies.vcf.bgz"
        this_vcf_filters_bgz="~{d}{this_basename}.reversed.withfilters.selfRef.homoplasmies.vcf.bgz"
        this_vcf_shifted_bgz="~{d}{this_basename}.reversed.selfRef.shifted.homoplasmies.vcf.bgz"


        bgzip -c "~{d}{this_mt_vcf}" > input_mt_vcf.gz && tabix input_mt_vcf.gz
        bcftools view -Oz -i 'FORMAT/AF>0.95' input_mt_vcf.gz > "~{d}{this_filt_vcf}.gz"
        tabix "~{d}{this_filt_vcf}.gz"
        bgzip -cd "~{d}{this_filt_vcf}.gz" > "~{d}{this_filt_vcf}"


        Rscript --vanilla ~{check_hom_overlap_script} "~{d}{this_sample}" "~{d}{this_filt_vcf}" "~{d}{this_nuc_vcf}" FALSE
        bgzip -c "overlapping_variants_to_remove.vcf" > mt_overlaps_rm.vcf.gz && tabix mt_overlaps_rm.vcf.gz
        bcftools isec "~{d}{this_filt_vcf}.gz" mt_overlaps_rm.vcf.gz -p output_isec
        export private_to_rem=$(cat ./output_isec/0001.vcf | grep ^chr | wc -l | sed 's/^ *//g')
        if [ $private_to_rem -ne 0 ]; then
          echo "ERROR: There should not be any variants private to the VCF for removal."
          exit 1;
        fi
        export shared_by_both=$(cat ./output_isec/0002.vcf | grep ^chr | wc -l | sed 's/^ *//g')
        export for_removal=$(cat overlapping_variants_to_remove.vcf | grep ^chr | wc -l | sed 's/^ *//g')
        if [ $shared_by_both -ne $for_removal ]; then
          echo "ERROR: The number of variants shared by both rejected VCF and MT VCF should be the number of variants for removal."
          exit 1;
        fi
        export only_undupe=$(cat ./output_isec/0000.vcf | grep ^chr | wc -l | sed 's/^ *//g')
        export original_nrow=$(cat "~{d}{this_filt_vcf}" | grep ^chr | wc -l | sed 's/^ *//g')
        if [ "$((original_nrow-for_removal))" -ne $only_undupe ]; then
          echo "ERROR: New VCF should be smaller than old VCF by the exact number of records for removal."
          exit 1;
        fi
        Rscript --vanilla ~{check_hom_overlap_script} "temp" "output_isec/0000.vcf" "~{d}{this_nuc_vcf}" TRUE
        rm "~{d}{this_filt_vcf}.gz" "~{d}{this_filt_vcf}" "~{d}{this_filt_vcf}.gz.tbi" mt_overlaps_rm.vcf.gz mt_overlaps_rm.vcf.gz.tbi overlapping_variants_to_remove.vcf
        cp "output_isec/0000.vcf" "~{d}{this_filt_vcf}"
        bgzip -c "~{d}{this_filt_vcf}" > "~{d}{this_filt_vcf}".gz && tabix "~{d}{this_filt_vcf}".gz
        rm -rf output_isec


        bcftools consensus -f mt_fasta_renamed.fasta -o mt_fasta_lifted.fasta -c "~{d}{this_mt_chain}" "~{d}{this_filt_vcf}.gz"
        Rscript --vanilla ~{fa_renaming_script} mt_fasta_lifted.fasta internal_mt.interval_list TRUE "~{d}{this_mt_fasta}" TRUE
        java -Xmx1000m -jar /usr/gitc/picard.jar \
          CreateSequenceDictionary \
          REFERENCE="~{d}{this_mt_fasta}" \
          OUTPUT="~{d}{this_basename}.dict"
        samtools faidx "~{d}{this_mt_fasta}"

        Rscript --vanilla ~{variant_bounds_script} "~{d}{this_nuc_vcf}" internal_nuc.interval_list rejected_nuc.vcf
        bgzip -c "~{d}{this_nuc_vcf}" > input_nuc_vcf.gz && tabix input_nuc_vcf.gz
        bgzip -c rejected_nuc.vcf > rejected_nuc.vcf.gz && tabix rejected_nuc.vcf.gz
        bcftools isec -p intersected_vcfs -Ov input_nuc_vcf.gz rejected_nuc.vcf.gz
        export private_to_rej=$(cat ./intersected_vcfs/0001.vcf | grep ^chr | wc -l | sed 's/^ *//g')
        if [ $private_to_rej -ne 0 ]; then
          echo "ERROR: There should not be any variants private to the rejected VCF."
          exit 1;
        fi
        cat ./intersected_vcfs/0002.vcf | grep ^chr | wc -l | sed 's/^ *//g' > "~{d}{this_sample}.nuc.removed.txt"

        java -jar /usr/gitc/picard.jar ExtractSequences \
          INTERVAL_LIST=internal_nuc.interval_list \
          R=~{ref_fasta} \
          O=nuc_fasta.fasta
        Rscript --vanilla ~{fa_renaming_script} nuc_fasta.fasta internal_nuc.interval_list FALSE nuc_fasta_renamed.fasta FALSE
        samtools faidx nuc_fasta_renamed.fasta
        bgzip -c ./intersected_vcfs/0000.vcf > filtered_nuc_vcf.gz && tabix filtered_nuc_vcf.gz
        bcftools consensus -f nuc_fasta_renamed.fasta -o nuc_fasta_lifted.fasta -c "~{d}{this_nuc_chain}" filtered_nuc_vcf.gz
        Rscript --vanilla ~{fa_renaming_script} nuc_fasta_lifted.fasta internal_nuc.interval_list TRUE "~{d}{this_nuc_only_fasta}" FALSE
        cat "~{d}{this_nuc_only_fasta}" "~{d}{this_mt_fasta}" > "~{d}{this_nuc_mt_fasta}"
        #/usr/gitc/bwa index ~{d}{this_nuc_mt_fasta}
        java -Xmx1000m -jar /usr/gitc/picard.jar \
          CreateSequenceDictionary \
          REFERENCE="~{d}{this_nuc_mt_fasta}" \
          OUTPUT="~{d}{this_basename_both}.dict"
        samtools faidx "~{d}{this_nuc_mt_fasta}"

        java -Xmx1000m -jar /usr/gitc/picard.jar \
          CreateSequenceDictionary \
          REFERENCE="~{d}{this_nuc_only_fasta}" \
          OUTPUT="~{d}{this_basename_nuc}.dict"
        samtools faidx "~{d}{this_nuc_only_fasta}"

        echo "Now shifting the reference..."
        R --vanilla <<CODE
          full_fasta <- readLines("~{d}{this_mt_fasta}")
          topline <- full_fasta[1]
          linelen <- nchar(full_fasta[2])
          n_shift <- ~{n_shift}
          other_lines <- paste0(full_fasta[2:length(full_fasta)],collapse='')
          other_lines_shifted <- paste0(substr(other_lines, n_shift+1, nchar(other_lines)), substr(other_lines, 1, n_shift))
          len_chr <- nchar(other_lines_shifted)
          shifted_split <- substring(other_lines_shifted, seq(1, len_chr, linelen), unique(c(seq(linelen, len_chr, linelen), len_chr)))
          final_data <- c(topline, shifted_split)
          writeLines(final_data, "~{d}{this_fasta_shifted}")

          total_len <- nchar(other_lines)
          sec1 <- paste(c('chain',9999,'chrM',total_len,'+', 0,total_len-n_shift, 'chrM', total_len, '+', n_shift, total_len, 1),collapse=' ')
          sec2 <- paste(c('chain',9999,'chrM',total_len,'+', total_len-n_shift,total_len, 'chrM', total_len, '+', 0, n_shift, 2),collapse=' ')
          writeLines(c(sec1, total_len-n_shift, '', sec2, n_shift, ''), "~{d}{this_chain_shifted}")

          sec1_f <- paste(c('chain',9999,'chrM',total_len,'+', n_shift,total_len, 'chrM', total_len, '+', 0, total_len-n_shift, 1),collapse=' ')
          sec2_f <- paste(c('chain',9999,'chrM',total_len,'+', 0, n_shift, 'chrM', total_len, '+', total_len-n_shift, total_len, 2),collapse=' ')
          writeLines(c(sec1_f, total_len-n_shift, '', sec2_f, n_shift, ''), "~{d}{this_chain_fwd_shifted}")
        CODE

        cat "~{d}{this_fasta_shifted}" "~{d}{this_nuc_only_fasta}" > "~{d}{this_fasta_cat_shifted}"

        #/usr/gitc/bwa index ~{d}{this_fasta_cat_shifted}
        java -Xmx1000m -jar /usr/gitc/picard.jar \
          CreateSequenceDictionary \
          REFERENCE="~{d}{this_fasta_cat_shifted}" \
          OUTPUT="~{d}{this_shifted_basename_append}.dict"
        samtools faidx "~{d}{this_fasta_cat_shifted}"

        java -Xmx1000m -jar /usr/gitc/picard.jar \
          CreateSequenceDictionary \
          REFERENCE="~{d}{this_fasta_shifted}" \
          OUTPUT="~{d}{this_shifted_basename}.dict"
        samtools faidx "~{d}{this_fasta_shifted}"

        echo "Now adjusting the MT interval list..."
        R --vanilla <<CODE
          intervals <- readLines("~{mt_interval_list}")
          intervals <- intervals[grep('^@', intervals, invert=T)]
          interval_names <- sapply(strsplit(intervals, '\t'),function(x)x[5])
          new_header <- readLines("~{d}{this_basename}.dict")
          lens <- sapply(interval_names, function(x)as.numeric(gsub('LN:','',strsplit(new_header[grep(paste0('^@SQ\tSN:',x), new_header)[1]], '\t')[[1]][3])))
          if(any(is.na(lens))) stop('ERROR: Some NUMT intervals were not found in the mt_andNuc sequence dictionary.')
          new_intervals <- c(new_header, paste(interval_names, 1, lens, '+', interval_names, sep='\t'))
          writeLines(new_intervals, "~{d}{this_mt_intervals}")
        CODE

        echo "Now shifting the noncontrol region..."
        java -jar /usr/gitc/picard.jar LiftOverIntervalList \
          I="~{non_control_region_interval_list}" \
          O="~{d}{this_nonctrl_interval}" \
          SD="~{d}{this_basename}.dict" \
          CHAIN="~{d}{this_mt_chain}"

        echo "Now shifting the control region..."
        R --vanilla <<CODE
          full_intervals <- readLines("~{d}{this_nonctrl_interval}")
          correct_dict <- readLines("~{d}{this_shifted_basename}.dict")
          n_shift <- ~{n_shift}
          if (length(full_intervals) > 3) {
            stop('ERROR: there should be only 3 lines in interval list.')
          }
          if ((length(grep('^@', full_intervals)) != 2) | (length(grep('^chrM', full_intervals)) != 1)) {
            stop('ERROR: there should be 2 comment lines and one interval line on chrM')
          }
          out_intervals <- c(full_intervals[1], correct_dict[2])
          split_line2 <- strsplit(out_intervals[2],'\t')[[1]]
          this_len <- as.numeric(gsub('^LN:','',split_line2[grep('^LN', split_line2)[1]]))

          interval <- strsplit(full_intervals[3],'\t')[[1]]
          new_end <- as.numeric(interval[2]) - 1 + this_len - n_shift
          new_start <- as.numeric(interval[3]) + 1 - n_shift
          interval[2:3] <- c(as.character(new_start), as.character(new_end))
          out_intervals <- c(out_intervals, paste0(interval, collapse='\t'))

          writeLines(out_intervals, "~{d}{this_ctrl_interval}")
        CODE

        echo "Now making force call variants..."
        python3.7 <<CODE
        import hail as hl

        def fai_to_len(fai):
            with open(fai) as f:
                line = f.readline()
            return int(line.split('\t')[1])

        def check_vcf_integrity(mt):
            # check that locus, alleles are the two key fields
            if sorted(list(mt.row_key)) != ['alleles', 'locus']:
                raise ValueError('VCFs must always be keyed by locus, alleles.')

            # check that all sites are bi-allelic
            if mt.aggregate_rows(~hl.agg.all(hl.len(mt.alleles) == 2)):
                raise ValueError('This function only supports biallelic sites (run SplitMultiAllelics!)')

            # check that there is no missingness in locus
            if mt.aggregate_rows(~hl.agg.all(hl.is_defined(mt.locus))):
                raise ValueError('ERROR: locus must always be defined, both before and after Liftover. This should be a reversible operation, thus finding missing loci after reverse liftover is very concerning.')

            # check that there is no missingness in alleles
            if mt.aggregate_rows(~hl.agg.all(hl.all(hl.map(hl.is_defined, mt.alleles)))):
                raise ValueError('ERROR: alleles should always be defined.')

        def apply_conversion(mt, liftover_target, skip_flip=False):
            mt_lifted = mt.annotate_rows(new_locus = hl.liftover(mt.locus, liftover_target))
            if skip_flip:
                mt_lifted = mt_lifted.annotate_rows(allele_flip = mt_lifted.alleles)
            else:
                mt_lifted = mt_lifted.annotate_rows(allele_flip = hl.reversed(mt_lifted.alleles))
            mt_lifted = mt_lifted.key_rows_by().rename({'locus':'locus_orig', 'alleles':'alleles_orig'}).rename({'new_locus':'locus', 'allele_flip': 'alleles'}).key_rows_by('locus','alleles')
            mt_lifted = mt_lifted.drop('locus_orig', 'alleles_orig')

            return mt_lifted

        target = hl.ReferenceGenome("target_self", ['chrM'], {'chrM':fai_to_len("~{d}{this_mt_fasta}.fai")}, mt_contigs=['chrM'])
        source = hl.ReferenceGenome('mtGRCh38', ['chrM'], {'chrM':fai_to_len("~{mt_ref_fasta_index}")}, mt_contigs=['chrM'])
        target.add_sequence("~{d}{this_mt_fasta}", "~{d}{this_mt_fasta}.fai")
        source.add_sequence("~{mt_ref_fasta}", "~{mt_ref_fasta_index}")
        source.add_liftover("~{d}{this_mt_chain}", "target_self")
        shifted_target = hl.ReferenceGenome("target_self_shifted", ['chrM'], {'chrM':fai_to_len("~{d}{this_fasta_shifted}.fai")}, mt_contigs=['chrM'])
        shifted_target.add_sequence("~{d}{this_fasta_shifted}", "~{d}{this_fasta_shifted}.fai")
        target.add_liftover("~{d}{this_chain_fwd_shifted}", shifted_target)

        mt_new = hl.import_vcf("~{d}{this_filt_vcf}", reference_genome='mtGRCh38').select_entries()
        check_vcf_integrity(mt_new)

        mt_new_1 = mt_new.select_rows()
        mt_lifted_target = apply_conversion(mt_new_1, "target_self")
        check_vcf_integrity(mt_lifted_target)
        hl.export_vcf(mt_lifted_target, "~{d}{this_vcf_bgz}", tabix=True)

        mt_new_2 = mt_new.select_rows('filters')
        this_metadata = hl.get_vcf_metadata("~{d}{this_filt_vcf}")
        this_metadata = {'filter': this_metadata['filter']}
        mt_lifted_target_withfilter = apply_conversion(mt_new_2, "target_self")
        check_vcf_integrity(mt_lifted_target_withfilter)
        hl.export_vcf(mt_lifted_target_withfilter, "~{d}{this_vcf_filters_bgz}", tabix=True, metadata=this_metadata)

        mt_lifted_shifted_target = apply_conversion(mt_lifted_target, "target_self_shifted", skip_flip=True)
        check_vcf_integrity(mt_lifted_shifted_target)
        hl.export_vcf(mt_lifted_shifted_target, "~{d}{this_vcf_shifted_bgz}", tabix=True)
        CODE
    >>>

    output {
        File self_fasta = "out/~{sample_name}~{suffix}.fasta"
        File self_fasta_index = "out/~{sample_name}~{suffix}.fasta.fai"
        File self_dict = "out/~{sample_name}~{suffix}.dict"

        File grch38_to_self_chain = "out/reference_to_~{sample_name}.chain"
        File filtered_vcf_ref_coord = "out/~{sample_name}.mtref.homoplasmies.vcf" # in reference coordinates

        File self_cat_fasta = "out/~{sample_name}andNuc~{suffix}.fasta"
        File self_cat_fasta_index = "out/~{sample_name}andNuc~{suffix}.fasta.fai"
        File self_cat_dict = "out/~{sample_name}andNuc~{suffix}.dict"

        File self_reference_nuc_fasta = "out/~{sample_name}NucOnly~{suffix}.fasta"
        File self_reference_nuc_fasta_index = "out/~{sample_name}NucOnly~{suffix}.fasta.fai"
        File self_reference_nuc_fasta_dict = "out/~{sample_name}NucOnly~{suffix}.dict"
        File grch38_to_self_nuc_chain = "out/reference_to_~{sample_name}NucOnly.chain"

        Int nuc_variants_dropped = read_int("out/~{sample_name}.nuc.removed.txt")
        Int mtdna_consensus_overlaps = read_int("out/~{sample_name}.mtdna_consensus_overlaps.txt")
        Int nuc_consensus_overlaps = read_int("out/~{sample_name}.nucdna_consensus_overlaps.txt")

        File self_shifted_fasta = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.fasta"
        File self_shifted_fasta_index = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.fasta.fai"
        File self_shifted_dict = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.dict"

        File shift_forward_chain = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.shift_fwd_~{n_shift}_bases.chain"
        File shift_back_chain = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.shift_back_~{n_shift}_bases.chain"

        File self_shifted_cat_fasta = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.cat.fasta"
        File self_shifted_cat_fasta_index = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.cat.fasta.fai"
        File self_shifted_cat_dict = "out/~{sample_name}~{suffix}.shifted_by_~{n_shift}_bases.cat.dict"

        File lifted_mt_intervals = "out/~{mt_intervals_basename}.~{sample_name}~{selfref_intervals_suffix}"
        File lifted_noncontrol_intervals = "out/~{noncontrol_basename}.~{sample_name}~{selfref_intervals_suffix}"
        File lifted_control_intervals = "out/control_region_shifted.chrM.~{sample_name}~{selfref_intervals_suffix}"

        File reversed_hom_vcf = "out/~{sample_name}~{suffix}.reversed.selfRef.homoplasmies.vcf.bgz"
        File reversed_hom_vcf_idx = "out/~{sample_name}~{suffix}.reversed.selfRef.homoplasmies.vcf.bgz.tbi"
        File reversed_hom_filters_vcf = "out/~{sample_name}~{suffix}.reversed.withfilters.selfRef.homoplasmies.vcf.bgz"
        File reversed_hom_filters_vcf_idx = "out/~{sample_name}~{suffix}.reversed.withfilters.selfRef.homoplasmies.vcf.bgz.tbi"
        File reversed_hom_vcf_shifted = "out/~{sample_name}~{suffix}.reversed.selfRef.shifted.homoplasmies.vcf.bgz"
        File reversed_hom_vcf_shifted_idx = "out/~{sample_name}~{suffix}.reversed.selfRef.shifted.homoplasmies.vcf.bgz.tbi"
    }

    runtime {
        docker: genomes_cloud_docker
        memory: "2 GB"
        disks: "local-disk 25 HDD"
        preemptible: select_first([preemptible_tries, 5])
    }
}

task MongoChainSwapLiftoverBed {
    input {
        File source_chain
        String input_target_name

        String input_source_name
        File input_bed
        File input_bed_index
        Int? preemptible_tries
        String ucsc_docker
    }

    Int disk_size = ceil(size(source_chain, "GB") + size(input_bed, "GB") * 2.5) * 4
    String this_bed_basename = basename(input_bed, '.bed')
    String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

    command <<<
        set -e

        mkdir out

        this_chain="~{source_chain}"
        this_target_sample="~{input_target_name}"
        this_target_path=out/"~{input_target_name}"
        this_base_path=out/"~{this_bed_basename}"

        chainSwap "~{d}{this_chain}" "~{d}{this_target_path}_to_~{input_source_name}.chain"
        liftOver ~{input_bed} "~{d}{this_chain}" "~{d}{this_base_path}.~{d}{this_target_sample}.liftedOver.bed" "~{d}{this_base_path}.~{d}{this_target_sample}.liftedOverRejects.unmapped"
        $igvtools index "~{d}{this_base_path}.~{d}{this_target_sample}.liftedOver.bed"
    >>>

    output {
        File chain = "out/~{input_target_name}_to_~{input_source_name}.chain"
        File transformed_bed = "out/~{this_bed_basename}.~{input_target_name}.liftedOver.bed"
        File transformed_bed_index = "out/~{this_bed_basename}.~{input_target_name}.liftedOver.bed.idx"
        File rejected = "out/~{this_bed_basename}.~{input_target_name}.liftedOverRejects.unmapped"
    }

    runtime {
        disks: "local-disk " + disk_size + " HDD"
        memory: "1 GB"
        docker: ucsc_docker
        preemptible: select_first([preemptible_tries, 5])
    }
}

task MongoHC {
    input {
        File ref_fasta
        File ref_fai
        File ref_dict

        File input_bam
        File input_bai

        String sample_name
        String suffix = ""

        Int max_reads_per_alignment_start = 75
        String? hc_extra_args
        Boolean make_bamout = false

        File? nuc_interval_list
        File? force_call_vcf
        File? force_call_vcf_index

        Boolean compress
        String gatk_version
        File? gatk_override
        String? gatk_docker_override
        Float? contamination

        Int hc_dp_lower_bound

        Int mem
        Int? preemptible_tries
        Int? n_cpu
    }

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB") + size(ref_dict, "GB")
    Int disk_size = ceil((size(input_bam, "GB") * 2) + ref_size) + 22

    String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        mkdir out
        this_sample=out/"~{sample_name}"
        this_basename="~{d}{this_sample}""~{suffix}"
        bamoutfile="~{d}{this_basename}.bamout.bam"
        # Cromwell doesn't like optional task outputs, so we have to touch this file
        touch "~{d}{bamoutfile}"

        if [[ ~{make_bamout} == 'true' ]]; then bamoutstr="--bam-output ~{d}{this_basename}.bamout.bam"; else bamoutstr=""; fi

        gatk --java-options "-Xmx~{command_mem}m" HaplotypeCaller \
          -R ~{ref_fasta} \
          -I ~{input_bam} \
          ~{"-L " + nuc_interval_list} \
          -O "~{d}{this_basename}.raw.vcf" \
          ~{hc_extra_args} \
          -contamination ~{default="0" contamination} \
          ~{"--genotype-filtered-alleles --alleles " + force_call_vcf} \
          --max-reads-per-alignment-start ~{max_reads_per_alignment_start} \
          --max-mnp-distance 0 \
          --annotation StrandBiasBySample \
          -G StandardAnnotation -G StandardHCAnnotation \
          -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 ~{d}{bamoutstr}

        echo "Now applying hard filters..."
        gatk --java-options "-Xmx~{command_mem}m" SelectVariants -V "~{d}{this_basename}.raw.vcf" -select-type SNP -O snps.vcf
        gatk --java-options "-Xmx~{command_mem}m" VariantFiltration -V snps.vcf \
          -R ~{ref_fasta} \
          -O snps_filtered.vcf \
          -filter "QD < 2.0" --filter-name "QD2" \
          -filter "QUAL < 30.0" --filter-name "QUAL30" \
          -filter "SOR > 3.0" --filter-name "SOR3" \
          -filter "FS > 60.0" --filter-name "FS60" \
          -filter "MQ < 40.0" --filter-name "MQ40" \
          -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
          -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
          --genotype-filter-expression "isHet == 1" --genotype-filter-name "isHetFilt" \
          --genotype-filter-expression "isHomRef == 1" --genotype-filter-name "isHomRefFilt" \
          ~{'--genotype-filter-expression "DP < ' + hc_dp_lower_bound + '" --genotype-filter-name "genoDP' + hc_dp_lower_bound + '"'}

        gatk --java-options "-Xmx~{command_mem}m" SelectVariants -V "~{d}{this_basename}.raw.vcf" -select-type INDEL -O indels.vcf
        gatk --java-options "-Xmx~{command_mem}m" VariantFiltration -V indels.vcf \
          -R ~{ref_fasta} \
          -O indels_filtered.vcf \
          -filter "QD < 2.0" --filter-name "QD2" \
          -filter "QUAL < 30.0" --filter-name "QUAL30" \
          -filter "FS > 200.0" --filter-name "FS200" \
          -filter "SOR > 10.0" --filter-name "SOR10" \
          -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
          --genotype-filter-expression "isHet == 1" --genotype-filter-name "isHetFilt" \
          --genotype-filter-expression "isHomRef == 1" --genotype-filter-name "isHomRefFilt" \
          ~{'--genotype-filter-expression "DP < ' + hc_dp_lower_bound + '" --genotype-filter-name "genoDP' + hc_dp_lower_bound + '"'}

        gatk --java-options "-Xmx~{command_mem}m" MergeVcfs -I snps_filtered.vcf -I indels_filtered.vcf -O "~{d}{this_basename}.vcf"

        echo "Now filtering VCF..."
        gatk --java-options "-Xmx~{command_mem}m" SelectVariants \
          -V "~{d}{this_basename}.vcf" \
          --exclude-filtered \
          --set-filtered-gt-to-nocall \
          --exclude-non-variants \
          -O "~{d}{this_basename}.pass.vcf"

        gatk --java-options "-Xmx~{command_mem}m" CountVariants -V $this_basename.pass.vcf | tail -n1 > "~{d}{this_basename}.passvars.txt"

        echo "Now splitting multi-allelics..."
        gatk --java-options "-Xmx~{command_mem}m" LeftAlignAndTrimVariants \
          -R ~{ref_fasta} \
          -V "~{d}{this_basename}.pass.vcf" \
          -O "~{d}{this_basename}.pass.split.vcf" \
          --split-multi-allelics \
          --dont-trim-alleles \
          --keep-original-ac \
          --create-output-variant-index
    >>>

    runtime {
        docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_tries, 5])
        cpu: select_first([n_cpu,1])
    }
    output {
        File raw_vcf = "out/~{sample_name}~{suffix}.raw.vcf"
        File raw_vcf_idx = "out/~{sample_name}~{suffix}.raw.vcf.idx"
        File output_bamOut = "out/~{sample_name}~{suffix}.bamout.bam"
        File filtered_vcf = "out/~{sample_name}~{suffix}.vcf"
        File filtered_vcf_idx = "out/~{sample_name}~{suffix}.vcf.idx"
        File full_pass_vcf = "out/~{sample_name}~{suffix}.pass.vcf"
        File full_pass_vcf_index = "out/~{sample_name}~{suffix}.pass.vcf.idx"
        Int post_filt_vars = read_int('out/~{sample_name}~{suffix}.passvars.txt')
        File split_vcf = "out/~{sample_name}~{suffix}.pass.split.vcf"
        File split_vcf_index = "out/~{sample_name}~{suffix}.pass.split.vcf.idx"
    }
}

task MongoNucM2 {
    input {
        File ref_fasta
        File ref_fai
        File ref_dict
        File input_bam
        File input_bai

        String sample_name
        String suffix = ""

        Int max_reads_per_alignment_start = 75
        String? m2_extra_args
        Boolean make_bamout = false
        Boolean compress

        File? mt_interval_list

        Float? vaf_cutoff
        String? m2_extra_filtering_args
        Int max_alt_allele_count
        Float? vaf_filter_threshold
        Float? f_score_beta
        Float? verifyBamID
        File? blacklisted_sites
        File? blacklisted_sites_index

        # runtime
        File? gatk_override
        String gatk_version
        String? gatk_docker_override
        Int mem
        Int? preemptible_tries
        Int? n_cpu
    }

    Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB")
    Int disk_size = ceil(size(input_bam, "GB")*2 + ref_size) + 20
    Float defval = 0.0

    # # hc_contamination will be None if hasContamination is not defined (I think) OR contamination_major not defined OR contamination_minor not defined
    # String hasContamination_2 = select_first([hasContamination,"NOT FOUND"])
    # Float? hc_contamination = if run_contamination && hasContamination_2 == "YES" then (if contamination_major == 0.0 then contamination_minor else 1.0 - contamination_major) else 0.0
    # Float hc_contamination_2 = select_first([hc_contamination, 0.0])
    # Float? max_contamination = if defined(verifyBamID) then (if verifyBamID > hc_contamination_2 then verifyBamID else hc_contamination_2) else hc_contamination_2

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

    meta {
        description: "Mutect2 for calling Snps and Indels"
    }
    parameter_meta {
        input_bam: "Aligned Bam"
    }

    #if [[ {defined(verifyBamID)} == 'true' ]]; then
    #  arr_contamination=('~{sep="' '" [select_first([verifyBamID, defval])]}')
    #else
    #  arr_contamination=$(printf '0.0 %.0s' {1..~{length([sample_name])}})
    #fi

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
        echo "Extra arguments for mutect2: ""~{m2_extra_args}""$cust_interval"

        mkdir out
        this_sample=out/"~{sample_name}"
        this_contamination="~{select_first([verifyBamID, defval])}"
        this_bam="~{input_bam}"
        this_basename="~{d}{this_sample}~{suffix}"
        bamoutfile="~{d}{this_basename}.bamout.bam"
        touch "~{d}{bamoutfile}"
        if [[ ~{make_bamout} == 'true' ]]; then bamoutstr="--bam-output ~{d}{bamoutfile}"; else bamoutstr=""; fi

        gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
          -R ~{ref_fasta} \
          -I "~{d}{this_bam}" \
          ~{"-L " + mt_interval_list} \
          -O "~{d}{this_basename}.raw.vcf" \
          ~{m2_extra_args} \
          ~{"--minimum-allele-fraction " + vaf_filter_threshold} \
          --annotation StrandBiasBySample \
          --max-reads-per-alignment-start ~{max_reads_per_alignment_start} \
          --max-mnp-distance 0 ~{d}{bamoutstr}

        echo "Now filtering..."
        gatk --java-options "-Xmx~{command_mem}m" FilterMutectCalls -V "~{d}{this_basename}.raw.vcf" \
          -R ~{ref_fasta} \
          -O filtered.vcf \
          --stats "~{d}{this_basename}.raw.vcf.stats" \
          ~{m2_extra_filtering_args} \
          --max-alt-allele-count ~{max_alt_allele_count} \
          ~{"--min-allele-fraction " + vaf_filter_threshold} \
          ~{"--f-score-beta " + f_score_beta} \
          --contamination-estimate "~{d}{this_contamination}"

        ~{"gatk IndexFeatureFile -I " + blacklisted_sites}

        gatk --java-options "-Xmx~{command_mem}m" VariantFiltration -V filtered.vcf \
          -O "~{d}{this_basename}.vcf" \
          --apply-allele-specific-filters \
          ~{"--mask-name 'blacklisted_site' --mask " + blacklisted_sites}

        echo "Now filtering VCF..."
        gatk --java-options "-Xmx~{command_mem}m" SelectVariants \
          -V "~{d}{this_basename}.vcf" \
          --exclude-filtered \
          -O "~{d}{this_basename}.pass.vcf"

        gatk CountVariants -V "~{d}{this_basename}.pass.vcf" | tail -n1 > "~{d}{this_basename}.passvars.txt"

        echo "Now splitting..."
        gatk --java-options "-Xmx~{command_mem}m" LeftAlignAndTrimVariants \
          -R ~{ref_fasta} \
          -V "~{d}{this_basename}.pass.vcf" \
          -O "~{d}{this_basename}.pass.split.vcf" \
          --split-multi-allelics \
          --dont-trim-alleles \
          --keep-original-ac \
          --create-output-variant-index
    >>>
    runtime {
        docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_tries, 5])
        cpu: select_first([n_cpu,2])
    }
    output {
        File raw_vcf = "out/~{sample_name}~{suffix}.raw.vcf"
        File raw_vcf_idx = "out/~{sample_name}~{suffix}.raw.vcf.idx"
        File stats = "out/~{sample_name}~{suffix}.raw.vcf.stats"
        File output_bamOut = "out/~{sample_name}~{suffix}.bamout.bam"

        File filtered_vcf = "out/~{sample_name}~{suffix}.vcf"
        File filtered_vcf_idx = "out/~{sample_name}~{suffix}.vcf.idx"

        File full_pass_vcf = "out/~{sample_name}~{suffix}.pass.vcf"
        File full_pass_vcf_index = "out/~{sample_name}~{suffix}.pass.vcf.idx"
        Int post_filt_vars = read_int('out/~{sample_name}~{suffix}.passvars.txt')

        File split_vcf = "out/~{sample_name}~{suffix}.pass.split.vcf"
        File split_vcf_index = "out/~{sample_name}~{suffix}.pass.split.vcf.idx"
    }
}

task MongoRunM2InitialFilterSplit {
    # intended for use prior to contamination (e.g., assumes mtDNA, no contamination)
    input {
        String sample_name
        File input_bam
        File input_bai
        Float? verifyBamID
        String suffix

        File ref_fasta
        File ref_fai
        File ref_dict
        Int max_reads_per_alignment_start = 75
        String? m2_extra_args
        Boolean make_bamout = false
        Boolean compress

        File? mt_interval_list

        Float? vaf_cutoff
        String? m2_extra_filtering_args
        Int max_alt_allele_count
        Float? vaf_filter_threshold
        Float? f_score_beta

        File? blacklisted_sites
        File? blacklisted_sites_index

        # runtime
        String? gatk_docker_override
        File? gatk_override
        String gatk_version
        Int mem
        Int? preemptible_tries
        Int? n_cpu
    }

    Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB")
    Int disk_size = (ceil(size(input_bam, "GB") + ref_size) * 2) + 20
    Float defval = 0.0

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

    meta {
        description: "Mutect2 for calling Snps and Indels + initial filter + split"
    }
    parameter_meta {
        input_bam: "Aligned Bam"
    }

    #if [[ {defined(verifyBamID)} == 'true' ]]; then
    #  arr_contamination=('~{sep="' '" [select_first([verifyBamID, defval])]}')
    #else
    #  arr_contamination=$(printf '0.0 %.0s' {1..~{length([sample_name])}})
    #fi

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
        echo "Extra arguments for mutect2: ""~{m2_extra_args}""$cust_interval"

        # We need to create these files regardless, even if they stay empty
        mkdir out
        this_sample=out/"~{sample_name}"
        this_contamination="~{select_first([verifyBamID, defval])}"
        this_basename="~{d}{this_sample}~{suffix}"
        bamoutfile="~{d}{this_basename}.bamout.bam"
        touch "~{d}{bamoutfile}"
        if [[ ~{make_bamout} == 'true' ]]; then bamoutstr="--bam-output ~{d}{bamoutfile}"; else bamoutstr=""; fi

        gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
          -R ~{ref_fasta} \
          -I ~{input_bam} \
          ~{"-L " + mt_interval_list} \
          -O "~{d}{this_basename}.raw.vcf" \
          ~{m2_extra_args} \
          --annotation StrandBiasBySample \
          --read-filter MateOnSameContigOrNoMappedMateReadFilter \
          --read-filter MateUnmappedAndUnmappedReadFilter \
          --mitochondria-mode \
          --max-reads-per-alignment-start ~{max_reads_per_alignment_start} \
          --max-mnp-distance 0 ~{d}{bamoutstr}

        echo "Now running filtering..."
        gatk --java-options "-Xmx~{command_mem}m" FilterMutectCalls -V "~{d}{this_basename}.raw.vcf" \
          -R ~{ref_fasta} \
          -O filtered.vcf \
          --stats "~{d}{this_basename}.raw.vcf.stats" \
          ~{m2_extra_filtering_args} \
          --max-alt-allele-count ~{max_alt_allele_count} \
          --mitochondria-mode \
          ~{"--min-allele-fraction " + vaf_filter_threshold} \
          ~{"--f-score-beta " + f_score_beta} \
          --contamination-estimate "~{d}{this_contamination}"

        ~{"gatk IndexFeatureFile -I " + blacklisted_sites}

        gatk --java-options "-Xmx~{command_mem}m" VariantFiltration -V filtered.vcf \
          -O "~{d}{this_basename}.filtered.vcf" \
          --apply-allele-specific-filters \
          ~{"--mask-name 'blacklisted_site' --mask " + blacklisted_sites}

        echo "Now splitting multi-allelics..."
        gatk --java-options "-Xmx~{command_mem}m" LeftAlignAndTrimVariants \
          -R ~{ref_fasta} \
          -V "~{d}{this_basename}.filtered.vcf" \
          -O split.vcf \
          --split-multi-allelics \
          --dont-trim-alleles \
          --keep-original-ac

        gatk --java-options "-Xmx~{command_mem}m" SelectVariants \
          -V split.vcf \
          -O "~{d}{this_basename}.splitAndPassOnly.vcf" \
          --exclude-filtered
    >>>
    runtime {
        docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_tries, 5])
        cpu: select_first([n_cpu,2])
    }
    output {
        File raw_vcf = "out/~{sample_name}~{suffix}.raw.vcf"
        File raw_vcf_idx = "out/~{sample_name}~{suffix}.raw.vcf.idx"
        File stats = "out/~{sample_name}~{suffix}.raw.vcf.stats"
        File output_bamOut = "out/~{sample_name}~{suffix}.bamout.bam"

        File filtered_vcf = "out/~{sample_name}~{suffix}.filtered.vcf"
        File filtered_vcf_idx = "out/~{sample_name}~{suffix}.filtered.vcf.idx"

        File vcf_for_haplochecker = "out/~{sample_name}~{suffix}.splitAndPassOnly.vcf"
    }
}

task MongoM2FilterContaminationSplit {
    input {
        File raw_vcf
        File raw_vcf_index
        File raw_vcf_stats
        String sample_name
        String hasContamination
        Float contamination_major
        Float contamination_minor
        Float? verifyBamID

        Boolean run_contamination
        File ref_fasta
        File ref_fai
        File ref_dict

        Boolean compress
        Float? vaf_cutoff
        String suffix

        String? m2_extra_filtering_args
        Int max_alt_allele_count
        Float? vaf_filter_threshold
        Float? f_score_beta

        File? blacklisted_sites
        File? blacklisted_sites_index

        File? gatk_override
        String? gatk_docker_override
        String gatk_version

        # runtime
        Int? preemptible_tries
    }

    Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB")
    Int disk_size = ceil(size(raw_vcf, "GB") + ref_size) + 20
    Float defval = 0.0
    String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

    # hc_contamination will be None if hasContamination is not defined (I think) OR contamination_major not defined OR contamination_minor not defined
    #Float hc_contamination = if hasContamination == "YES" then (if contamination_major == 0.0 then contamination_minor else 1.0 - contamination_major) else 0.0
    #Float? max_contamination = if defined(verifyBamID) then (if verifyBamID > hc_contamination then verifyBamID else hc_contamination) else hc_contamination

    meta {
        description: "Mutect2 Filtering for calling Snps and Indels"
    }
    parameter_meta {
        vaf_filter_threshold: "Hard cutoff for minimum allele fraction. All sites with VAF less than this cutoff will be filtered."
        f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
    }

    #if [[ {defined(verifyBamID)} == 'true' ]]; then
    #  arr_verifybam_contamination=('~{sep="' '" [select_first([verifyBamID, defval])]}')
    #else
    #  arr_verifybam_contamination=$(printf '0.0 %.0s' {1..~{length([sample_name])}})
    #fi

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        mkdir out

        this_sample=out/"~{sample_name}"
        this_raw_vcf="~{raw_vcf}"
        this_raw_stats="~{raw_vcf_stats}"
        this_has_contam="~{hasContamination}"
        this_verifybam="~{select_first([verifyBamID, defval])}"
        this_contam_major="~{contamination_major}"
        this_contam_minor="~{contamination_minor}"

        this_basename="~{d}{this_sample}~{suffix}"
        bamoutfile="~{d}{this_basename}.bamout.bam"
        touch "~{d}{bamoutfile}"

        if [[ "~{d}{this_has_contam}" == 'YES' ]]; then
          if (( $(echo "~{d}{this_contam_major} == 0.0"|bc -l) )); then
            this_hc_contamination="~{d}{this_contam_minor}"
          else
            this_hc_contamination=$( bc <<< "1-~{d}{this_contam_major}" )
          fi
        else
          this_hc_contamination=0.0
        fi

        echo "~{d}{this_hc_contamination}" > "~{d}{this_basename}.hc_contam.txt"

        if (( $(echo "~{d}{this_verifybam} > ~{d}{this_hc_contamination}"|bc -l) )); then
          this_max_contamination="~{d}{this_verifybam}"
        else
          this_max_contamination="~{d}{this_hc_contamination}"
        fi

        echo "VerifyBam: ~{d}{this_verifybam}; HC: ~{d}{this_hc_contamination}; Max: ~{d}{this_max_contamination}"

        gatk --java-options "-Xmx2500m" FilterMutectCalls \
          -V "~{d}{this_raw_vcf}" \
          -R ~{ref_fasta} \
          -O filtered.vcf \
          --stats "~{d}{this_raw_stats}" \
          ~{m2_extra_filtering_args} \
          --max-alt-allele-count ~{max_alt_allele_count} \
          --mitochondria-mode \
          ~{"--min-allele-fraction " + vaf_filter_threshold} \
          ~{"--f-score-beta " + f_score_beta} \
          --contamination-estimate "~{d}{this_max_contamination}"

        ~{"gatk IndexFeatureFile -I " + blacklisted_sites}

        gatk --java-options "-Xmx2500m" VariantFiltration \
          -V filtered.vcf \
          -O "~{d}{this_basename}.vcf" \
          --apply-allele-specific-filters \
          ~{"--mask-name 'blacklisted_site' --mask " + blacklisted_sites}

        echo "Now splitting multi-allelics..."
        gatk --java-options "-Xmx2500m" LeftAlignAndTrimVariants \
          -R ~{ref_fasta} \
          -V "~{d}{this_basename}.vcf" \
          -O "~{d}{this_basename}.split.vcf" \
          --split-multi-allelics \
          --dont-trim-alleles \
          --keep-original-ac \
          --create-output-variant-index

    >>>
    runtime {
        docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
        memory: "4 MB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_tries, 5])
        cpu: 2
    }
    output {
        File filtered_vcf = "out/~{sample_name}~{suffix}.vcf"
        File filtered_vcf_idx = "out/~{sample_name}~{suffix}.vcf.idx"
        File split_vcf = "out/~{sample_name}~{suffix}.split.vcf"
        File split_vcf_index = "out/~{sample_name}~{suffix}.split.vcf.idx"
        Float contamination = read_float("out/~{sample_name}~{suffix}.hc_contam.txt") # now an optional output, producing UNDEFINED if not computed
    }
}

task MongoAlignToMtRegShiftedAndMetrics {
    input {
        File input_bam
        String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -Y $bash_ref_fasta"
        String sample_base_name
        String suffix

        File mt
        File mt_index
        File mt_dict

        File mt_cat
        File mt_cat_index
        File mt_cat_dict

        File mt_shifted
        File mt_shifted_index
        File mt_shifted_dict

        File mt_shifted_cat
        File mt_shifted_cat_index
        File mt_shifted_cat_dict

        File mt_interval_list

        String? read_name_regex

        Int? read_length
        Int? coverage_cap
        Int? n_cpu

        Int? preemptible_tries
    }

    Int this_cpu = select_first([n_cpu, 2])
    String this_bwa_commandline = bwa_commandline + " -t " + this_cpu
    Float ref_size = size(mt_cat, "GB") + size(mt_cat_index, "GB")
    Float shifted_ref_size = size(mt_shifted_cat, "GB") + size(mt_shifted_cat_index, "GB")
    Int disk_size = ceil(size(input_bam, "GB") * 4 + ref_size * 4 + shifted_ref_size * 4) + 20

    Int read_length_for_optimization = select_first([read_length, 151])
    String this_output_bam_basename = basename(input_bam, '.bam') + ".remap" + suffix

    String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

    meta {
        description: "Aligns with BWA and MergeBamAlignment, then Marks Duplicates. Outputs a coordinate sorted bam."
    }
    parameter_meta {
        input_bam: "Unmapped bam"
    }
    command <<<

        export BWAVERSION=$(/usr/gitc/bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')

        set -o pipefail
        set -e

        mkdir out

        this_sample=out/"~{sample_base_name}~{suffix}"
        this_bam="~{input_bam}"
        this_output_bam_basename=out/"$(basename ~{d}{this_bam} .bam).remap~{suffix}"
        this_mt_intervals="~{mt_interval_list}"
        this_mt_cat_fasta="~{mt_cat}"
        this_mt_fasta="~{mt}"
        this_mt_shifted_cat_fasta="~{mt_shifted_cat}"
        this_mt_shifted_fasta="~{mt_shifted}"

        # set the bash variable needed for the command-line
        /usr/gitc/bwa index "~{d}{this_mt_cat_fasta}"
        bash_ref_fasta="~{d}{this_mt_cat_fasta}"
        java -Xms5000m -jar /usr/gitc/picard.jar \
          SamToFastq \
          INPUT="~{d}{this_bam}" \
          FASTQ=/dev/stdout \
          INTERLEAVE=true \
          NON_PF=true | \
        /usr/gitc/~{this_bwa_commandline} /dev/stdin - 2> >(tee "~{d}{this_output_bam_basename}.bwa.stderr.log" >&2) | \
        java -Xms5000m -jar /usr/gitc/picard.jar \
          MergeBamAlignment \
          VALIDATION_STRINGENCY=SILENT \
          EXPECTED_ORIENTATIONS=FR \
          ATTRIBUTES_TO_RETAIN=X0 \
          ATTRIBUTES_TO_REMOVE=NM \
          ATTRIBUTES_TO_REMOVE=MD \
          ALIGNED_BAM=/dev/stdin \
          UNMAPPED_BAM="~{d}{this_bam}" \
          OUTPUT=mba.bam \
          REFERENCE_SEQUENCE="~{d}{this_mt_cat_fasta}" \
          PAIRED_RUN=true \
          SORT_ORDER="unsorted" \
          IS_BISULFITE_SEQUENCE=false \
          ALIGNED_READS_ONLY=false \
          CLIP_ADAPTERS=false \
          MAX_RECORDS_IN_RAM=2000000 \
          ADD_MATE_CIGAR=true \
          MAX_INSERTIONS_OR_DELETIONS=-1 \
          PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
          PROGRAM_RECORD_ID="bwamem" \
          PROGRAM_GROUP_VERSION=$BWAVERSION \
          PROGRAM_GROUP_COMMAND_LINE="~{this_bwa_commandline}" \
          PROGRAM_GROUP_NAME="bwamem" \
          UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
          ALIGNER_PROPER_PAIR_FLAGS=true \
          UNMAP_CONTAMINANT_READS=true \
          ADD_PG_TAG_TO_READS=false

        java -Xms5000m -jar /usr/gitc/picard.jar \
          MarkDuplicates \
          INPUT=mba.bam \
          OUTPUT=md.bam \
          METRICS_FILE="~{d}{this_output_bam_basename}.metrics" \
          VALIDATION_STRINGENCY=SILENT \
          ~{"READ_NAME_REGEX=" + read_name_regex} \
          OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
          ASSUME_SORT_ORDER="queryname" \
          CLEAR_DT="false" \
          ADD_PG_TAG_TO_READS=false

        java -Xms5000m -jar /usr/gitc/picard.jar \
          SortSam \
          INPUT=md.bam \
          OUTPUT="~{d}{this_output_bam_basename}_pre_mt_filt.bam" \
          SORT_ORDER="coordinate" \
          CREATE_INDEX=true \
          MAX_RECORDS_IN_RAM=300000

        # now we have to subset to mito and update sequence dictionary
        java -Xms5000m -jar /usr/gitc/picard.jar \
          ReorderSam \
          I="~{d}{this_output_bam_basename}_pre_mt_filt.bam" \
          O="~{d}{this_output_bam_basename}.bam" \
          REFERENCE="~{d}{this_mt_fasta}" \
          ALLOW_INCOMPLETE_DICT_CONCORDANCE=true \
          CREATE_INDEX=true

        echo "Now starting on shifted..."
        # set the bash variable needed for the command-line
        /usr/gitc/bwa index "~{d}{this_mt_shifted_cat_fasta}"
        bash_ref_fasta="~{d}{this_mt_shifted_cat_fasta}"
        java -Xms5000m -jar /usr/gitc/picard.jar \
          SamToFastq \
          INPUT="~{d}{this_bam}" \
          FASTQ=/dev/stdout \
          INTERLEAVE=true \
          NON_PF=true | \
        /usr/gitc/~{this_bwa_commandline} /dev/stdin - 2> >(tee "~{d}{this_output_bam_basename}.shifted.bwa.stderr.log" >&2) | \
        java -Xms5000m -jar /usr/gitc/picard.jar \
          MergeBamAlignment \
          VALIDATION_STRINGENCY=SILENT \
          EXPECTED_ORIENTATIONS=FR \
          ATTRIBUTES_TO_RETAIN=X0 \
          ATTRIBUTES_TO_REMOVE=NM \
          ATTRIBUTES_TO_REMOVE=MD \
          ALIGNED_BAM=/dev/stdin \
          UNMAPPED_BAM="~{d}{this_bam}" \
          OUTPUT=mba.shifted.bam \
          REFERENCE_SEQUENCE="~{d}{this_mt_shifted_cat_fasta}" \
          PAIRED_RUN=true \
          SORT_ORDER="unsorted" \
          IS_BISULFITE_SEQUENCE=false \
          ALIGNED_READS_ONLY=false \
          CLIP_ADAPTERS=false \
          MAX_RECORDS_IN_RAM=2000000 \
          ADD_MATE_CIGAR=true \
          MAX_INSERTIONS_OR_DELETIONS=-1 \
          PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
          PROGRAM_RECORD_ID="bwamem" \
          PROGRAM_GROUP_VERSION=$BWAVERSION \
          PROGRAM_GROUP_COMMAND_LINE="~{this_bwa_commandline}" \
          PROGRAM_GROUP_NAME="bwamem" \
          UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
          ALIGNER_PROPER_PAIR_FLAGS=true \
          UNMAP_CONTAMINANT_READS=true \
          ADD_PG_TAG_TO_READS=false

        java -Xms5000m -jar /usr/gitc/picard.jar \
          MarkDuplicates \
          INPUT=mba.shifted.bam \
          OUTPUT=md.shifted.bam \
          METRICS_FILE="~{d}{this_output_bam_basename}.shifted.metrics" \
          VALIDATION_STRINGENCY=SILENT \
          ~{"READ_NAME_REGEX=" + read_name_regex} \
          OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
          ASSUME_SORT_ORDER="queryname" \
          CLEAR_DT="false" \
          ADD_PG_TAG_TO_READS=false

        java -Xms5000m -jar /usr/gitc/picard.jar \
          SortSam \
          INPUT=md.shifted.bam \
          OUTPUT="~{d}{this_output_bam_basename}.shifted_pre_mt_filt.bam" \
          SORT_ORDER="coordinate" \
          CREATE_INDEX=true \
          MAX_RECORDS_IN_RAM=300000

        # now we have to subset to mito and update sequence dictionary
        java -Xms5000m -jar /usr/gitc/picard.jar \
          ReorderSam \
          I="~{d}{this_output_bam_basename}.shifted_pre_mt_filt.bam" \
          O="~{d}{this_output_bam_basename}.shifted.bam" \
          REFERENCE="~{d}{this_mt_shifted_fasta}" \
          ALLOW_INCOMPLETE_DICT_CONCORDANCE=true \
          CREATE_INDEX=true

        echo "Now collecting wgs metrics..."
        java -Xms5000m -jar /usr/gitc/picard.jar \
          CollectWgsMetrics \
          INPUT="~{d}{this_output_bam_basename}.bam" \
          INTERVALS="~{d}{this_mt_intervals}" \
          VALIDATION_STRINGENCY=SILENT \
          REFERENCE_SEQUENCE="~{d}{this_mt_fasta}" \
          OUTPUT="~{d}{this_sample}_r2_wgs_metrics.txt" \
          USE_FAST_ALGORITHM=true \
          READ_LENGTH=~{read_length_for_optimization} \
          ~{"COVERAGE_CAP=" + coverage_cap} \
          INCLUDE_BQ_HISTOGRAM=true \
          THEORETICAL_SENSITIVITY_OUTPUT="~{d}{this_sample}_r2_wgs_theoretical_sensitivity.txt"

        R --vanilla <<CODE
          df = read.table("~{d}{this_sample}_r2_wgs_metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
          write.table(floor(df[,"MEAN_COVERAGE"]), "~{d}{this_sample}_r2_mean_coverage.txt", quote=F, col.names=F, row.names=F)
          write.table(df[,"MEDIAN_COVERAGE"], "~{d}{this_sample}_r2_median_coverage.txt", quote=F, col.names=F, row.names=F)
        CODE
    >>>
    runtime {
        preemptible: select_first([preemptible_tries, 5])
        memory: "6 GB"
        cpu: this_cpu
        disks: "local-disk " + disk_size + " HDD"
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
    }
    output {
        File mt_aligned_bam = "out/~{this_output_bam_basename}.bam"
        File mt_aligned_bai = "out/~{this_output_bam_basename}.bai"
        File nuc_and_mt_aligned_bam = "out/~{this_output_bam_basename}_pre_mt_filt.bam"
        File nuc_and_mt_aligned_bai = "out/~{this_output_bam_basename}_pre_mt_filt.bai"
        File bwa_stderr_log = "out/~{this_output_bam_basename}.bwa.stderr.log"
        File duplicate_metrics = "out/~{this_output_bam_basename}.metrics"

        File wgs_metrics = "out/~{sample_base_name}~{suffix}_r2_wgs_metrics.txt"
        File theoretical_sensitivity = "out/~{sample_base_name}~{suffix}_r2_wgs_theoretical_sensitivity.txt"
        Int mean_coverage = read_int("out/~{sample_base_name}~{suffix}_r2_mean_coverage.txt")
        Float median_coverage = read_float("out/~{sample_base_name}~{suffix}_r2_median_coverage.txt")

        File shifted_mt_aligned_bam = "out/~{this_output_bam_basename}.shifted.bam"
        File shifted_mt_aligned_bai = "out/~{this_output_bam_basename}.shifted.bai"
        File nuc_and_shifted_mt_aligned_bam = "out/~{this_output_bam_basename}.shifted_pre_mt_filt.bam"
        File nuc_and_shifted_mt_aligned_bai = "out/~{this_output_bam_basename}.shifted_pre_mt_filt.bai"
        File shifted_bwa_stderr_log = "out/~{this_output_bam_basename}.shifted.bwa.stderr.log"
        File shifted_duplicate_metrics = "out/~{this_output_bam_basename}.shifted.metrics"
    }
}

task MongoCallMtAndShifted {
    input {
        String sample_base_name
        File mt_self
        File mt_self_index
        File mt_self_dict
        File input_bam
        File input_bai
        File mt_interval_list
        File force_call_vcf
        File force_call_vcf_idx

        File shifted_mt_self
        File shifted_mt_self_index
        File shifted_mt_self_dict
        File shifted_input_bam
        File shifted_input_bai
        File shifted_mt_interval_list
        File shifted_force_call_vcf
        File shifted_force_call_vcf_idx

        String? m2_extra_args
        String? shifted_m2_extra_args

        Int max_reads_per_alignment_start = 75
        Boolean make_bamout = false
        Boolean compress
        String suffix

        # runtime
        File? gatk_override
        String gatk_version
        String? gatk_docker_override
        Int mem
        Int? preemptible_tries
        Int? n_cpu
    }

    Float ref_size = size(mt_self, "GB") + size(mt_self_index, "GB") + size(shifted_mt_self, "GB") + size(shifted_mt_self_index, "GB")
    Int disk_size = ceil(size(input_bam, "GB") + size(shifted_input_bam, "GB") + ref_size) + 20

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

    meta {
        description: "Mutect2 for calling Snps and Indels; runs on both shifted and nonshifted"
    }
    parameter_meta {
        input_bam: "Aligned Bam"
    }
    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
        echo "Extra arguments for mutect2: ""~{m2_extra_args}""$cust_interval"

        mkdir out

        this_sample=out/"~{sample_base_name}~{suffix}"

        this_bam="~{input_bam}"
        this_noncontrol="~{mt_interval_list}"
        this_force_vcf="~{force_call_vcf}"
        this_self_fasta="~{mt_self}"
        this_shifted_bam="~{shifted_input_bam}"
        this_control="~{shifted_mt_interval_list}"
        this_shifted_force_vcf="~{shifted_force_call_vcf}"
        this_self_shifted_fasta="~{shifted_mt_self}"

        touch "~{d}{this_sample}.bamout.bam"
        touch "~{d}{this_sample}.shifted.bamout.bam"

        if [[ ~{make_bamout} == 'true' ]]; then bamoutstr="--bam-output ~{d}{this_sample}.bamout.bam"; else bamoutstr=""; fi
        if [[ ~{make_bamout} == 'true' ]]; then shiftedbamoutstr="--bam-output ~{d}{this_sample}.shifted.bamout.bam"; else shiftedbamoutstr=""; fi

        echo "Obtaining force calls for specified VCF: ~{d}{this_force_vcf}"

        # Fix for DNANexus weirdness
        gatk IndexFeatureFile -I "~{d}{this_force_vcf}"

        gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
          -R "~{d}{this_self_fasta}" \
          -I "~{d}{this_bam}" \
          -L "~{d}{this_noncontrol}" \
          -O "~{d}{this_sample}.raw.vcf" \
          --genotype-filtered-alleles \
          --alleles "~{d}{this_force_vcf}" \
          ~{m2_extra_args} \
          --annotation StrandBiasBySample \
          --mitochondria-mode \
          --read-filter MateOnSameContigOrNoMappedMateReadFilter \
          --read-filter MateUnmappedAndUnmappedReadFilter \
          --max-reads-per-alignment-start ~{max_reads_per_alignment_start} \
          --max-mnp-distance 0 ~{d}{bamoutstr}

        echo "Obtaining force calls for specified VCF: ~{d}{this_shifted_force_vcf}"

        # Fix for DNANexus weirdness
        gatk IndexFeatureFile -I "~{d}{this_shifted_force_vcf}"

        gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
          -R "~{d}{this_self_shifted_fasta}" \
          -I "~{d}{this_shifted_bam}" \
          -L "~{d}{this_control}" \
          -O "~{d}{this_sample}.shifted.raw.vcf" \
          --genotype-filtered-alleles \
          --alleles "~{d}{this_shifted_force_vcf}" \
          ~{shifted_m2_extra_args} \
          --annotation StrandBiasBySample \
          --mitochondria-mode \
          --read-filter MateOnSameContigOrNoMappedMateReadFilter \
          --read-filter MateUnmappedAndUnmappedReadFilter \
          --max-reads-per-alignment-start ~{max_reads_per_alignment_start} \
          --max-mnp-distance 0 ~{d}{shiftedbamoutstr}

    >>>
    runtime {
        docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_tries, 5])
        cpu: select_first([n_cpu,2])
    }
    output {
        File raw_vcf = "out/~{sample_base_name}~{suffix}.raw.vcf"
        File raw_vcf_idx = "out/~{sample_base_name}~{suffix}.raw.vcf.idx"
        File stats = "out/~{sample_base_name}~{suffix}.raw.vcf.stats"
        File output_bamOut = "out/~{sample_base_name}~{suffix}.bamout.bam"

        File shifted_raw_vcf = "out/~{sample_base_name}~{suffix}.shifted.raw.vcf"
        File shifted_raw_vcf_idx = "out/~{sample_base_name}~{suffix}.shifted.raw.vcf.idx"
        File shifted_stats = "out/~{sample_base_name}~{suffix}.shifted.raw.vcf.stats"
        File shifted_output_bamOut = "out/~{sample_base_name}~{suffix}.shifted.bamout.bam"
    }
}

task MongoLiftoverCombineMergeFilterContamSplit {
    input {
        String sample_base_name
        String suffix

        File mt_self
        File mt_self_index
        File mt_self_dict
        File shifted_vcf
        File shifted_vcf_idx
        File non_shifted_vcf
        File non_shifted_vcf_idx
        File shifted_stats
        File non_shifted_stats
        File shift_back_chain

        String hasContamination
        Float contamination_major
        Float contamination_minor
        Float? verifyBamID
        File blacklisted_sites
        File blacklisted_sites_index

        Boolean compress
        Float? vaf_cutoff
        String? m2_extra_filtering_args
        Int max_alt_allele_count
        Float? vaf_filter_threshold
        Float? f_score_beta

        Int? preemptible_tries
        File? gatk_override
        String? gatk_docker_override
        String gatk_version
    }

    Float ref_size = size(mt_self, "GB") + size(mt_self_index, "GB")
    Int disk_size = ceil(size(shifted_vcf, "GB")*4 + ref_size) + 30
    Float defval = 0.0
    String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

    meta {
        description: "Lifts over shifted vcf of control region and combines it with the rest of the chrM calls and also merges stats."
    }
    parameter_meta {
        shifted_vcf: "VCF of control region on shifted reference"
        non_shifted_vcf: "VCF of the rest of chrM on original reference"
        mt_self: "Original (not shifted) chrM reference"
        shift_back_chain: "Chain file to lift over from shifted reference to original chrM"
    }

    #if [[ {defined(verifyBamID)} == 'true' ]]; then
    #  arr_verifybam_contamination=('~{sep="' '" [select_first([verifyBamID, defval])]}')
    #else
    #  arr_verifybam_contamination=$(printf '0.0 %.0s' {1..~{length([sample_base_name])}})
    #fi
    # Note that the perl below removes the ##contig=chrM style line from the file
    command<<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        mkdir out

        this_sample=out/"~{sample_base_name}~{suffix}"
        this_shifted_vcf="~{shifted_vcf}"
        this_nonshifted_vcf="~{non_shifted_vcf}"
        this_ref_fasta="~{mt_self}"
        this_shiftback_chain="~{shift_back_chain}"
        this_shifted_stats="~{shifted_stats}"
        this_nonshifted_stats="~{non_shifted_stats}"
        this_blacklisted="~{blacklisted_sites}"

        this_has_contam="~{hasContamination}"
        this_verifybam="~{select_first([verifyBamID, defval])}"
        this_contam_major="~{contamination_major}"
        this_contam_minor="~{contamination_minor}"

        if [[ "~{d}{this_has_contam}" == 'YES' ]]; then
          if (( $(echo "~{d}{this_contam_major} == 0.0"|bc -l) )); then
            this_hc_contamination="~{d}{this_contam_minor}"
          else
            this_hc_contamination=$( bc <<< "1-~{d}{this_contam_major}" )
          fi
        else
          this_hc_contamination=0.0
        fi

        echo "~{d}{this_hc_contamination}" > "~{d}{this_sample}.r2.hc_contam.txt"

        if (( $(echo "~{d}{this_verifybam} > ~{d}{this_hc_contamination}"|bc -l) )); then
          this_max_contamination="~{d}{this_verifybam}"
        else
          this_max_contamination="~{d}{this_hc_contamination}"
        fi

        echo "VerifyBam: ~{d}{this_verifybam}; HC: ~{d}{this_hc_contamination}; Max: ~{d}{this_max_contamination}"

        this_basename_raw="~{d}{this_sample}.raw"
        touch "~{d}{this_sample}.merge.bamout.bam"

        gatk LiftoverVcf \
          -I "~{d}{this_shifted_vcf}" \
          -O "~{d}{this_basename_raw}.shifted_back.vcf" \
          -R "~{d}{this_ref_fasta}" \
          -CHAIN "~{d}{this_shiftback_chain}" \
          -REJECT "~{d}{this_basename_raw}.rejected.vcf"

        cat "~{d}{this_nonshifted_vcf}" | perl -ne 'print unless /^##contig=<ID=(?!chrM)/' > chrMFiltVCF.vcf

        gatk MergeVcfs \
          -I "~{d}{this_basename_raw}.shifted_back.vcf" \
          -I chrMFiltVCF.vcf \
          -O "~{d}{this_basename_raw}.merged.vcf"

        cat "~{d}{this_basename_raw}.rejected.vcf" | grep ^chrM | wc -l > "~{d}{this_sample}.n_failed_vars.txt"

        gatk MergeMutectStats \
          --stats "~{d}{this_shifted_stats}" \
          --stats "~{d}{this_nonshifted_stats}" \
          -O "~{d}{this_basename_raw}.combined.stats"

        echo "Now filtering contamination..."
        gatk --java-options "-Xmx1000m" FilterMutectCalls -V "~{d}{this_basename_raw}.merged.vcf" \
          -R "~{d}{this_ref_fasta}" \
          -O filtered.vcf \
          --stats "~{d}{this_basename_raw}.combined.stats" \
          ~{m2_extra_filtering_args} \
          --max-alt-allele-count ~{max_alt_allele_count} \
          --mitochondria-mode \
          ~{"--min-allele-fraction " + vaf_filter_threshold} \
          ~{"--f-score-beta " + f_score_beta} \
          --contamination-estimate "~{d}{this_max_contamination}"

        gatk IndexFeatureFile -I "~{d}{this_blacklisted}"

        gatk --java-options "-Xmx1000m" VariantFiltration -V filtered.vcf \
          -O "~{d}{this_sample}.vcf" \
          --apply-allele-specific-filters \
          --mask-name 'blacklisted_site' \
          --mask "~{d}{this_blacklisted}"

        echo "Now splitting multi-allelics..."
        gatk --java-options "-Xmx1000m" LeftAlignAndTrimVariants \
          -R "~{d}{this_ref_fasta}" \
          -V "~{d}{this_sample}.vcf" \
          -O "~{d}{this_sample}.split.vcf" \
          --split-multi-allelics \
          --dont-trim-alleles \
          --keep-original-ac \
          --create-output-variant-index
    >>>
    output {
        File rejected_raw_vcf = "out/~{sample_base_name}~{suffix}.raw.rejected.vcf"
        File merged_raw_vcf = "out/~{sample_base_name}~{suffix}.raw.merged.vcf"
        File merged_raw_vcf_idx = "out/~{sample_base_name}~{suffix}.raw.merged.vcf.idx"
        Int number_failed = read_int("out/~{sample_base_name}~{suffix}.n_failed_vars.txt")
        File raw_stats = "out/~{sample_base_name}~{suffix}.raw.combined.stats"

        File filtered_vcf = "out/~{sample_base_name}~{suffix}.vcf"
        File filtered_vcf_idx = "out/~{sample_base_name}~{suffix}.vcf.idx"
        File split_vcf = "out/~{sample_base_name}~{suffix}.split.vcf"
        File split_vcf_idx = "out/~{sample_base_name}~{suffix}.split.vcf.idx"
        Float contamination = read_float('out/~{sample_base_name}~{suffix}.r2.hc_contam.txt') # now an optional output, producing UNDEFINED if not computed
    }
    runtime {
        docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:"+gatk_version])
        memory: "1200 MB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_tries, 5])
    }
}

task MongoLiftoverVCFAndGetCoverage {
    # A specialized routine to return the resultant VCF back to GRCh38
    input {
        String sample_name
        File original_filtered_vcf
        File new_self_ref_vcf
        File reversed_hom_ref_vcf

        File mt_self
        File mt_self_index
        File mt_self_dict

        File mt_self_shifted
        File mt_self_shifted_index
        File mt_self_shifted_dict

        File chain_self_to_ref
        File chain_ref_to_self

        File input_bam_regular_ref
        File input_bam_regular_ref_index
        File input_bam_shifted_ref
        File input_bam_shifted_ref_index
        File self_control_region_shifted_reference_interval_list
        File self_non_control_region_interval_list

        String self_suffix
        File HailLiftover
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        # runtime
        Int? n_cpu
        Int? preemptible_tries
        String genomes_cloud_docker
    }

    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(mt_self, "GB") + size(mt_self_shifted, "GB")
    Float bam_size = size(input_bam_regular_ref, 'GB') + size(input_bam_shifted_ref, 'GB')
    Int disk_size = ceil(bam_size) + ceil(size(original_filtered_vcf, "GB") + size(new_self_ref_vcf, "GB") + size(reversed_hom_ref_vcf, "GB") + ref_size) *2 + 20
    String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

    command <<<
        set -e

        mkdir out

        this_sample_name="~{sample_name}"
        this_self_ref_vcf="~{new_self_ref_vcf}"
        this_ref_filtered_vcf="~{original_filtered_vcf}"
        this_rev_hom_ref_vcf="~{reversed_hom_ref_vcf}"
        this_self_to_ref_chain="~{chain_self_to_ref}"
        this_ref_to_self_chain="~{chain_ref_to_self}"
        this_self_fasta="~{mt_self}"
        this_self_fai="~{mt_self_index}"
        this_self_shifted_fasta="~{mt_self_shifted}"
        this_self_control="~{self_control_region_shifted_reference_interval_list}"
        this_self_non_control="~{self_non_control_region_interval_list}"
        this_self_bam="~{input_bam_regular_ref}"
        this_self_shifted_bam="~{input_bam_shifted_ref}"

        this_sample=out/"~{d}{this_sample_name}"
        this_basename="~{d}{this_sample}~{self_suffix}.split"
        this_logging="~{d}{this_basename}_fix_liftover.log"

        bgzip -c "~{d}{this_self_ref_vcf}" > "~{d}{this_self_ref_vcf}.bgz" && tabix "~{d}{this_self_ref_vcf}.bgz"
        tabix "~{d}{this_rev_hom_ref_vcf}"
        bcftools isec -p intersected_vcfs -Ov "~{d}{this_self_ref_vcf}.bgz" "~{d}{this_rev_hom_ref_vcf}"

        # there should be no records private to reversed hom ref VCF
        export private_to_rev_hom_ref=$(cat ./intersected_vcfs/0001.vcf | grep ^chrM | wc -l | sed 's/^ *//g')
        if [ $private_to_rev_hom_ref -ne 0 ]; then
          echo "ERROR: There should not be any variants private to the reversed hom ref VCF."
          exit 1;
        fi

        java -jar /usr/gitc/picard.jar LiftoverVcf \
          I=./intersected_vcfs/0000.vcf \
          O="~{d}{this_basename}.selfToRef.pre.vcf" \
          R=~{ref_fasta} \
          CHAIN="~{d}{this_self_to_ref_chain}" \
          REJECT="~{d}{this_basename}.selfToRef.rejected.pre.vcf"

        java -jar /usr/gitc/picard.jar MergeVcfs \
          I="~{d}{this_basename}.selfToRef.rejected.pre.vcf" \
          I=./intersected_vcfs/0002.vcf \
          O="~{d}{this_basename}.selfToRef.rejected.vcf"

        sed -e 's/^chr//' "~{d}{this_ref_filtered_vcf}" \
          | awk '{OFS="\t"; if (!/^#/){print $1,$2-(length($4)>1 ? 0 : 1),$2-1+length($4),$4"/"$5,"+",length($4),length($5)}}' > filtered.bed
        sed -e 's/^chr//' "~{d}{this_basename}.selfToRef.pre.vcf" \
          | awk '{OFS="\t"; if (!/^#/){print $1,$2-(length($4)>1 ? 0 : 1),$2-1+length($4),$4"/"$5,"+",length($4),length($5)}}' > success.bed
        export n_ref_pass_thru=$(bedtools intersect -a filtered.bed -b success.bed | awk '{OFS="\t"; if (($7 == 1)) {print}}' | wc -l | sed 's/^ *//g')
        echo $n_ref_pass_thru > "~{d}{this_sample}_n_ref_pass_thru.txt"

        if [ $n_ref_pass_thru -ne 0 ]; then
          echo "ERROR: All variants changed in the self-reference, and all sites within self-reference insertions (excluding the first base), should have failed LiftoverVcf, which is not the case here."
          exit 1;
        fi

        export n_filtered=$(cat "~{d}{this_ref_filtered_vcf}" | grep ^chrM | wc -l | sed 's/^ *//g')
        export n_original=$(cat "~{d}{this_self_ref_vcf}" | grep ^chrM | wc -l | sed 's/^ *//g')
        export n_pass=$(cat "~{d}{this_basename}.selfToRef.pre.vcf" | grep ^chrM | wc -l | sed 's/^ *//g')
        export n_rejected=$(cat "~{d}{this_basename}.selfToRef.rejected.vcf" | grep ^chrM | wc -l | sed 's/^ *//g')
        export n_pass_rejected=$((n_pass + n_rejected))
        echo $n_pass > "~{d}{this_sample}_n_pass.txt"

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
        --vcf-file "~{d}{this_basename}.selfToRef.rejected.vcf" \
        --success-vcf-file "~{d}{this_basename}.selfToRef.pre.vcf" \
        --self-homoplasmies "~{d}{this_rev_hom_ref_vcf}" \
        --individual-name "~{d}{this_sample_name}" \
        --self-to-ref-chain "~{d}{this_self_to_ref_chain}" \
        --ref-to-self-chain "~{d}{this_ref_to_self_chain}" \
        --self-fasta "~{d}{this_self_fasta}" \
        --self-fai "~{d}{this_self_fai}" \
        --reference-fasta ~{ref_fasta} \
        --reference-fai ~{ref_fasta_index} \
        --output-prefix "~{d}{this_basename}.round2liftover" \
        --export-homoplasmic-deletions-coverage \
        --output-txt-for-wdl \
        --logging "~{d}{this_logging}"

        bgzip -cd "~{d}{this_basename}.round2liftover.rejected.vcf.bgz" > "~{d}{this_basename}.round2liftover.rejected.vcf"
        bgzip -cd "~{d}{this_basename}.round2liftover.fixed.vcf.bgz" > "~{d}{this_basename}.round2liftover.fixed.vcf"
        bgzip -cd "~{d}{this_basename}.round2liftover.updated_success.vcf.bgz" > "~{d}{this_basename}.round2liftover.updated_success.vcf"

        java -jar /usr/gitc/picard.jar MergeVcfs \
          I="~{d}{this_basename}.round2liftover.updated_success.vcf" \
          I="~{d}{this_basename}.round2liftover.fixed.vcf" \
          O="~{d}{this_basename}.selfToRef.final.vcf"

        export n_final_pass=$(cat "~{d}{this_basename}.selfToRef.final.vcf" | grep ^chrM | wc -l | sed 's/^ *//g')
        echo $n_final_pass > "~{d}{this_sample}_n_final_pass.txt"

        echo "Now producing coverage file..."
        java -jar /usr/gitc/picard.jar CollectHsMetrics \
          I="~{d}{this_self_bam}" \
          R="~{d}{this_self_fasta}" \
          PER_BASE_COVERAGE=non_control_region.tsv \
          O=non_control_region.metrics \
          TI="~{d}{this_self_non_control}" \
          BI="~{d}{this_self_non_control}" \
          COVMAX=20000 \
          SAMPLE_SIZE=1

        java -jar /usr/gitc/picard.jar CollectHsMetrics \
          I="~{d}{this_self_shifted_bam}" \
          R="~{d}{this_self_shifted_fasta}" \
          PER_BASE_COVERAGE=control_region_shifted.tsv \
          O=control_region_shifted.metrics \
          TI="~{d}{this_self_control}" \
          BI="~{d}{this_self_control}" \
          COVMAX=20000 \
          SAMPLE_SIZE=1

        R --vanilla <<CODE
          full_fasta <- readLines("~{d}{this_self_fasta}") # edited to account for variable reference sizes
          nlen <- nchar(paste0(full_fasta[2:length(full_fasta)],collapse=''))
          nshift <- 8000
          shift_back <- function(x) {
            if (x < (nlen-nshift+1)) {
              return(x + nshift)
            } else {
              return (x - (nlen-nshift))
          }
        }

        control_region_shifted = read.table("control_region_shifted.tsv", header=T)
        shifted_back = sapply(control_region_shifted[,"pos"], shift_back)
        control_region_shifted[,"pos"] = shifted_back

        beginning = subset(control_region_shifted, control_region_shifted[,'pos']<8000)
        end = subset(control_region_shifted, control_region_shifted[,'pos']>8000)

        non_control_region = read.table("non_control_region.tsv", header=T)
        combined_table = rbind(beginning, non_control_region, end)
        write.table(combined_table, "~{d}{this_basename}.per_base_coverage.tsv", row.names=F, col.names=T, quote=F, sep="\t")

        CODE

        echo "Now outputting a final table with integers..."
        paste -d "\t" "~{d}{this_basename}.round2liftover.all_int_outputs.txt" <(printf "n_liftover_changed_selfref_and_passed\n$(cat ~{d}{this_sample}_n_ref_pass_thru.txt)\n") <(printf "n_liftover_r1_pass\n$(cat ~{d}{this_sample}_n_pass.txt)\n") <(printf "n_liftover_r2_pass\n$(cat ~{d}{this_sample}_n_final_pass.txt)\n") > "~{d}{this_basename}.round2liftover.all_int_outputs.final.txt"
    >>>

    runtime {
        disks: "local-disk " + disk_size + " HDD"
        memory: "3000 MB"
        cpu: select_first([n_cpu, 2])
        docker: genomes_cloud_docker
        preemptible: select_first([preemptible_tries, 5])
    }

    output {
        File liftover_r1_rejected_vcf = "out/~{sample_name}~{self_suffix}.split.selfToRef.rejected.vcf"
        File liftover_r1_vcf = "out/~{sample_name}~{self_suffix}.split.selfToRef.pre.vcf"
        File liftover_r2_success_r1_vcf = "out/~{sample_name}~{self_suffix}.split.round2liftover.updated_success.vcf"
        File liftover_r2_rejected_vcf = "out/~{sample_name}~{self_suffix}.split.round2liftover.rejected.vcf"
        File liftover_r2_intermediate_vcf = "out/~{sample_name}~{self_suffix}.split.round2liftover.fixed.vcf"
        File liftover_r2_final_vcf = "out/~{sample_name}~{self_suffix}.split.selfToRef.final.vcf"
        File liftover_r2_log = "out/~{sample_name}~{self_suffix}.split_fix_liftover.log"
        File gap_coverage = "out/~{sample_name}~{self_suffix}.split.round2liftover.deletions_coverage.tsv"
        File self_coverage_table = "out/~{sample_name}~{self_suffix}.split.per_base_coverage.tsv"
        File liftoverStats = "out/~{sample_name}~{self_suffix}.split.round2liftover.all_int_outputs.final.txt"

        # stats
        Int n_liftover_changed_selfref_and_passed = read_int('out/~{sample_name}_n_ref_pass_thru.txt')
        Int n_liftover_r2_failed = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.round2_failed_sites.txt')
        Int n_liftover_r2_fixed = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.round2_fixed_sites.txt')
        Int n_liftover_r2_pass = read_int('out/~{sample_name}_n_final_pass.txt')
        Int n_liftover_r2_left_shift = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.left_alignment_of_indels.txt')
        Int n_liftover_r2_injected_from_success = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.round2_success_injected.txt')
        Int n_liftover_r2_ref_insertion_new_haplo = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.ref_insertion_new_haplos.txt')
        Int n_liftover_r2_failed_het_dele_span_insertion_boundary = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.het_deletions_span_insertions.txt')
        Int n_liftover_r2_failed_new_dupes_leftshift = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.new_dupes_left_shift_failed.txt')
        Int n_liftover_r2_het_ins_sharing_lhs_hom_dele = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.het_insertions_sharing_lhs_with_hom_ref_deletion.txt')
        Int n_liftover_r2_spanning_complex = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.het_dele_span_insert_repaired_with_complex_rework.txt')
        Int n_liftover_r2_spanningfixrhs_sharedlhs = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.heteroplasmic_deletions_sharing_lhs_with_homoplasmic_insertion_spanning_rhs.txt')
        Int n_liftover_r2_spanningfixlhs_upstream = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.heteroplasmic_deletions_spanning_only_lhs_of_homoplasmic_insertion.txt')
        Int n_liftover_r2_repaired_success = read_int('out/~{sample_name}~{self_suffix}.split.round2liftover.success_sites_flipped.txt')
    }
}

task MongoLiftoverSelfAndCollectOutputs {
    # run liftover provided a coverage txt file
    # use reference allele depth from the SELF-REFERENCE data after reversing the force-call
    # to impute the appropriate coverage for a deleted site. This REF allele is the
    # actual original sequence.
    input {
        String sample_name
        File self_ref_table
        File chain
        File homoplasmic_deletions_coverage
        File liftover_table

        Int mean_coverage
        Float median_coverage
        String major_haplogroup
        Float contamination
        Int nuc_variants_pass
        Int n_reads_unpaired_dropped
        Int nuc_variants_dropped
        Int mtdna_consensus_overlaps
        Int nuc_consensus_overlaps

        String ucsc_docker
        Int? preemptible_tries
    }

    Int disk_size = ceil(size(self_ref_table, "GB") * 4) + 20
    String d = "$" # a stupid trick to get ${} indexing in bash to work in Cromwell

    command <<<
        set -e

        mkdir out

        this_sample=out/"~{sample_name}"
        this_chain="~{chain}"
        this_dele_suppl="~{homoplasmic_deletions_coverage}"
        this_self_ref_tab="~{self_ref_table}"

        this_litover_stats="~{liftover_table}"
        this_mean_cov="~{mean_coverage}"
        this_median_cov="~{median_coverage}"
        this_maj_haplo="~{major_haplogroup}"
        this_contam="~{contamination}"
        this_nuc_var_pass="~{nuc_variants_pass}"
        this_n_reads_drop="~{n_reads_unpaired_dropped}"
        this_nuc_var_drop="~{nuc_variants_dropped}"
        this_mt_overlap="~{mtdna_consensus_overlaps}"
        this_nuc_overlap="~{nuc_consensus_overlaps}"

        tail -n +2 "~{d}{this_self_ref_tab}" | awk '{OFS="\t"} {print $1"\t"$2-1"\t"$2"\t"$4}' > per_base_coverage.bed
        liftOver per_base_coverage.bed "~{d}{this_chain}" "~{d}{this_sample}.liftedOver.bed" "~{d}{this_sample}.liftedOver.unmapped"

        R --vanilla <<CODE
          init_chain <- readLines("~{d}{this_chain}")
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

          chain <- read.csv("~{d}{this_chain}", skip=1, header=F, col.names=c('bases_in_block', 'ns', 'nt'), sep=sep, stringsAsFactors=F, fill=T)
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
          failure_bed <- read.csv("~{d}{this_sample}.liftedOver.unmapped", header=F, sep='\t', col.names=c('chr','start','end','cov'), comment.char='#')
          if (n_rows_different_keys_bed(expected_failed_liftover, failure_bed) > 0) {
            stop('ERROR: Rejected bed file is not identical to expectation.')
          }

          target_len <- as.numeric(strsplit(init_chain[1], '\\\\s+')[[1]][9])
          deletions <- chain[chain[['nt']] > 0,]
          expected_liftover_gaps <- start_end_to_bed(deletions[['coord_t']], deletions[['coord_t']]+deletions[['nt']])
          success_bed <- read.csv("~{d}{this_sample}.liftedOver.bed", header=F, sep='\t', col.names=c('chr','start','end','cov'), comment.char='#')
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
            fill_in_coverage <- read.csv("~{d}{this_dele_suppl}", sep='\t', stringsAsFactors=F)
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
          write.table(final_tsv[,c('chrom','pos','target','coverage')], "~{d}{this_sample}.appended.liftedOver.tsv", row.names=F, sep='\t', quote=F)
        CODE

        R --vanilla <<CODE
          dfLiftover = read.csv("~{d}{this_litover_stats}", sep='\t', stringsAsFactors=F)
          df = data.frame(mean_coverage = "~{d}{this_mean_cov}",
                          median_coverage = "~{d}{this_median_cov}",
                          major_haplogroup = "~{d}{this_maj_haplo}",
                          contamination = "~{d}{this_contam}",
                          nuc_variants_pass = "~{d}{this_nuc_var_pass}",
                          n_reads_unpaired_dropped = "~{d}{this_n_reads_drop}",
                          nuc_variants_dropped = "~{d}{this_nuc_var_drop}",
                          mtdna_consensus_overlaps = "~{d}{this_mt_overlap}",
                          nuc_consensus_overlaps = "~{d}{this_nuc_overlap}")
          write.table(cbind(dfLiftover, df), "~{d}{this_sample}_mtanalysis_diagnostic_statistics.tsv", row.names=F, sep='\t', quote=F)
        CODE
    >>>

    output {
        File reference_coverage = "out/~{sample_name}.appended.liftedOver.tsv"
        File intermediate_bed = "out/~{sample_name}.liftedOver.bed"
        File rejected = "out/~{sample_name}.liftedOver.unmapped"
        File table = "out/~{sample_name}_mtanalysis_diagnostic_statistics.tsv"
    }

    runtime {
        disks: "local-disk " + disk_size + " HDD"
        memory: "2 GB"
        docker: ucsc_docker
        preemptible: select_first([preemptible_tries, 5])
    }
}