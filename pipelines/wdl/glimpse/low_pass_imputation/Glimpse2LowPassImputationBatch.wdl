version 1.0

workflow Glimpse2LowPassImputationBatch {
    String pipeline_version = "1.1.0"

    input {
        Array[String] contigs
        Array[Array[String]] chromosome_groups # e.g., [["chr1","chr2","chr3"], ["chr4","chr5","chr6","chr7"], ...] must cover all contigs continuously in order
        
        String reference_panel_prefix

        File cram_manifest
        File fasta
        File fasta_index
        String output_basename

        Boolean impute_reference_only_variants = false
        Boolean call_indels = false

        Map[String, Array[String]] bcftools_shard_map

        Int? glimpse_phase_cpu_override

        String gatk_docker
        String glimpse_docker
    }

    Int defined_glimpse_phase_cpu_override = select_first([glimpse_phase_cpu_override, 4])

    call ParseCramManifest {
        input:
            cram_manifest = cram_manifest
    }

    # Step 1: Pre-calculate counts (Localizing whole-genome reference files here is removed for efficiency)
    scatter(contig_index in range(length(contigs))) {
        String current_contig = contigs[contig_index]
        Array[String] current_bcftools_shards = bcftools_shard_map[current_contig]
        Int current_shard_count = length(current_bcftools_shards)
    }

    # Step 2: Calculate global shard indices and perfectly subset the localization arrays per group
    call CalculateShardIndices {
        input:
            shard_counts = current_shard_count,
            chromosome_groups = chromosome_groups,
            reference_panel_prefix = reference_panel_prefix
    }

    # Step 3: Extract utilizing explicit Contig Groups with highly-optimized disk/CPU loops
    scatter(sample_idx in range(length(ParseCramManifest.crams))) {
        scatter(group_idx in range(length(chromosome_groups))) {
            
            call ExtractGenotypeLikelihoods {
                input:
                    cram = ParseCramManifest.crams[sample_idx],
                    cram_index = ParseCramManifest.cram_indices[sample_idx],
                    chromosome_group = chromosome_groups[group_idx],
                    contig_shard_starts = CalculateShardIndices.starts,
                    contigs = contigs,
                    flat_shards = flatten(current_bcftools_shards),
                    shard_counts = current_shard_count,
                    sites_vcfs = CalculateShardIndices.group_sites_vcfs[group_idx],
                    sites_vcf_indices = CalculateShardIndices.group_sites_vcf_indices[group_idx],
                    sites_tables = CalculateShardIndices.group_sites_tables[group_idx],
                    sites_table_indices = CalculateShardIndices.group_sites_table_indices[group_idx],
                    fasta = fasta,
                    fasta_index = fasta_index,
                    call_indels = call_indels
            }
        }
        
        # Flattens the scattered BCF outputs and indices into single sequential arrays for this sample
        Array[File] sample_flat_bcfs = flatten(ExtractGenotypeLikelihoods.flat_bcfs)
        Array[File] sample_flat_bcf_indices = flatten(ExtractGenotypeLikelihoods.flat_bcf_indices)
    }

    # WDL Native Transpose: Converts [Samples][Total_Shards] -> [Total_Shards][Samples] instantly
    Array[Array[File]] global_shard_x_sample = transpose(sample_flat_bcfs)
    Array[Array[File]] global_shard_x_sample_indices = transpose(sample_flat_bcf_indices)

    # Step 4: Index Arithmetic loop to process the FLAT merge natively
    scatter(contig_index in range(length(contigs))) {
        String current_contig_for_merge = contigs[contig_index]
        Array[Int] global_indices_for_contig = CalculateShardIndices.indices[contig_index]

        scatter (local_shard_idx in range(length(global_indices_for_contig))) {
            Int global_shard_idx = global_indices_for_contig[local_shard_idx]
            
            Array[File] shard_sample_bcfs = global_shard_x_sample[global_shard_idx]
            Array[File] shard_sample_bcf_indices = global_shard_x_sample_indices[global_shard_idx]
            
            # Execute standard bcftools merge natively using co-located files
            call BcftoolsMerge {
                input:
                    bcfs = shard_sample_bcfs,
                    bcf_indices = shard_sample_bcf_indices,
                    output_basename = output_basename + "." + current_contig_for_merge + ".shard_" + local_shard_idx
            }
        }

        call BcftoolsConcatNaive {
            input:
                bcfs = BcftoolsMerge.merged_bcf,
                output_basename = output_basename + "." + current_contig_for_merge + ".concat"
        }
    }

    output {
        Int total_samples = ParseCramManifest.total_samples
    }
}

task CalculateShardIndices {
    input {
        Array[Int] shard_counts
        Array[Array[String]] chromosome_groups
        String reference_panel_prefix
    }

    command <<<
        python3 <<EOF
        import json
        
        counts = [~{sep="," shard_counts}]
        indices = []
        starts = []
        curr = 0
        for c in counts:
            indices.append(list(range(curr, curr + c)))
            starts.append(curr)
            curr += c
            
        with open("indices.json", "w") as f:
            json.dump(indices, f)
        with open("starts.json", "w") as f:
            json.dump(starts, f)
            
        # Dynamically calculate exact paths to prevent massive array localization upstream
        prefix = "~{reference_panel_prefix}"
        with open("~{write_json(chromosome_groups)}", "r") as f:
            groups = json.load(f)
            
        vcfs, vcf_idxs, tbls, tbl_idxs = [], [], [], []
        for g in groups:
            chroms = [c.strip() for c in g]
            vcfs.append([f"{prefix}sites.{c}.vcf.gz" for c in chroms])
            vcf_idxs.append([f"{prefix}sites.{c}.vcf.gz.tbi" for c in chroms])
            tbls.append([f"{prefix}sites_table.{c}.gz" for c in chroms])
            tbl_idxs.append([f"{prefix}sites_table.{c}.gz.tbi" for c in chroms])
            
        with open("vcfs.json", "w") as f: json.dump(vcfs, f)
        with open("vcf_idxs.json", "w") as f: json.dump(vcf_idxs, f)
        with open("tbls.json", "w") as f: json.dump(tbls, f)
        with open("tbl_idxs.json", "w") as f: json.dump(tbl_idxs, f)
        EOF
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk 10 SSD"
        memory: "1 GiB"
        preemptible: 3
        noAddress: true
    }

    output {
        Array[Array[Int]] indices = read_json("indices.json")
        Array[Int] starts = read_json("starts.json")
        
        Array[Array[String]] group_sites_vcfs = read_json("vcfs.json")
        Array[Array[String]] group_sites_vcf_indices = read_json("vcf_idxs.json")
        Array[Array[String]] group_sites_tables = read_json("tbls.json")
        Array[Array[String]] group_sites_table_indices = read_json("tbl_idxs.json")
    }
}

task ParseCramManifest {
    input {
        File cram_manifest
    }

    command <<<
        cat <<EOF > script.py
        import json
        import pandas as pd
        import sys

        df = pd.read_csv("~{cram_manifest}", sep='\t')

        required_cols = ['cram_path', 'cram_index_path']
        missing_cols = [col for col in required_cols if col not in df.columns]

        if missing_cols:
            print(f"Missing required columns in the CRAM manifest: {', '.join(missing_cols)}.", file=sys.stderr)
            exit(1)

        crams = df['cram_path'].tolist()
        cram_indices = df['cram_index_path'].tolist()

        with open('crams.json', 'w') as json_file:
            json.dump(crams, json_file)
        with open('cram_indices.json', 'w') as json_file:
            json.dump(cram_indices, json_file)

        with open('total_samples.txt', 'w') as total_sample_file:
            total_sample_file.write(str(len(crams)))
        EOF
        python3 script.py
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        cpu: 1
        disks: "local-disk 10 SSD"
        memory: "1 GiB"
        preemptible: 3
        noAddress: true
    }

    output {
        Array[String] crams = read_json('crams.json')
        Array[String] cram_indices = read_json('cram_indices.json')
        Int total_samples = read_int("total_samples.txt")
    }
}

task ExtractGenotypeLikelihoods {
    input {
        File cram
        File cram_index
        
        Array[String] chromosome_group
        Array[Int] contig_shard_starts
        
        Array[String] contigs
        Array[String] flat_shards
        Array[Int] shard_counts
        
        Array[File] sites_vcfs
        Array[File] sites_vcf_indices
        Array[File] sites_tables
        Array[File] sites_table_indices
        
        File fasta
        File fasta_index
        Boolean call_indels = false

        Int seed = 12345
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 1
    }

    Int disk_size_gb = ceil(1.5 * size(cram, "GiB") + size(fasta, "GiB") + size(sites_vcfs, "GiB") + size(sites_tables, "GiB") + 10)

    command <<<
        set -xeuo pipefail

        ln -sf ~{cram} "~{basename(cram)}"
        ln -sf ~{cram_index} "~{basename(cram)}.crai"

        CONTIGS=(~{sep=" " contigs})
        VCFS=(~{sep=" " sites_vcfs})
        VCF_IDXS=(~{sep=" " sites_vcf_indices})
        TABLES=(~{sep=" " sites_tables})
        TABLE_IDXS=(~{sep=" " sites_table_indices})
        SHARD_COUNTS=(~{sep=" " shard_counts})
        SHARD_STARTS=(~{sep=" " contig_shard_starts})

        for i in "${!VCFS[@]}"; do
            ln -sf "${VCFS[$i]}" "vcf_${i}.vcf.gz"
            ln -sf "${VCF_IDXS[$i]}" "vcf_${i}.vcf.gz.tbi"
            ln -sf "${TABLES[$i]}" "tbl_${i}.gz"
            ln -sf "${TABLE_IDXS[$i]}" "tbl_${i}.gz.tbi"
        done

        ALL_SHARDS=()
        while IFS= read -r line || [ -n "$line" ]; do
            clean_line="${line//[$'\r\n']/}"
            if [ -n "$clean_line" ]; then
                ALL_SHARDS+=("$clean_line")
            fi
        done < "~{write_lines(flat_shards)}"

        GROUP_CONTIGS=(~{sep=" " chromosome_group})
        
        g_idx=0
        for raw_chrom in "${GROUP_CONTIGS[@]}"; do
            chrom="${raw_chrom// /}"
            
            i=-1
            for idx in "${!CONTIGS[@]}"; do
                if [[ "${CONTIGS[$idx]}" == "$chrom" ]]; then
                    i=$idx
                    break
                fi
            done
            
            if [[ $i -eq -1 ]]; then
                echo "Error: Contig $chrom not found in global contigs array!" >&2
                exit 1
            fi
            
            # Map correctly to localized subset arrays via g_idx
            vcf="vcf_${g_idx}.vcf.gz"
            tbl="tbl_${g_idx}.gz"
            count="${SHARD_COUNTS[$i]}"
            shard_idx="${SHARD_STARTS[$i]}"
            pad_contig=$(printf "%03d" $i)
            
            # 1. mpileup whole contig at once to prevent FASTA/VCF re-parsing
            bcftools mpileup --no-version -f ~{fasta} ~{if !call_indels then "-I " else ""} --seed ~{seed} -E -r "${chrom}" -T "${vcf}" -Ou ~{basename(cram)} | \
            bcftools call --no-version -Aim -C alleles -T "${tbl}" -Ou | \
            bcftools annotate --no-version -x 'INFO,^FORMAT/PL' -Ou | \
            bcftools norm --no-version -m -both -Ob -o "contig_${pad_contig}.bcf"
            
            bcftools index "contig_${pad_contig}.bcf"
            
            # 2. Slice master BCF into requested shards or rename if unsharded (1 total shard)
            if [[ $count -eq 1 ]]; then
                out_bcf="out_c${pad_contig}_s0000.bcf"
                
                mv "contig_${pad_contig}.bcf" "${out_bcf}"
                mv "contig_${pad_contig}.bcf.csi" "${out_bcf}.csi"
                
                ((++shard_idx))
            else
                for (( j=0; j<count; j++ )); do
                    shard="${ALL_SHARDS[$shard_idx]}"
                    pad_shard=$(printf "%04d" $j)
                    out_bcf="out_c${pad_contig}_s${pad_shard}.bcf"
                    
                    bcftools view -r "${shard}" -Ob -o "${out_bcf}" "contig_${pad_contig}.bcf"
                    bcftools index "${out_bcf}"
                    
                    ((++shard_idx))
                done
                
                rm "contig_${pad_contig}.bcf" "contig_${pad_contig}.bcf.csi"
            fi
            
            ((++g_idx))
        done
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        disks: "local-disk " + disk_size_gb + " SSD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
        noAddress: true
    }

    output {
        Array[File] flat_bcfs = glob("out_c*_s*.bcf")
        Array[File] flat_bcf_indices = glob("out_c*_s*.bcf.csi")
    }
}

task BcftoolsMerge {
    input {
        Array[File] bcfs
        Array[File] bcf_indices
        Int mem_gb = 6
        Int cpu = 1
        Int preemptible = 0
        Int max_retries = 3

        String output_basename
    }

    Int disk_size_gb = ceil(3 * size(bcfs, "GiB")) + 50

    command <<<
        set -euo pipefail
        
        BCFS=(~{sep=" " bcfs})
        IDXS=(~{sep=" " bcf_indices})
        
        > file_list.txt
        
        for i in "${!BCFS[@]}"; do
            # Prefix the symlink with the array index to prevent basename collisions
            safe_bcf="${i}_$(basename "${BCFS[$i]}")"
            safe_idx="${i}_$(basename "${IDXS[$i]}")"
            
            ln -sf "${BCFS[$i]}" "$safe_bcf"
            ln -sf "${IDXS[$i]}" "$safe_idx"
            
            # Sync filesystem timestamps to prevent htslib older-index warnings
            touch "$safe_idx"
            
            echo "$safe_bcf" >> file_list.txt
        done
        
        bcftools merge --no-version -O b -o ~{output_basename}.merged.bcf -l file_list.txt
        bcftools index ~{output_basename}.merged.bcf
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
        noAddress: true
    }

    output {
        File merged_bcf = "~{output_basename}.merged.bcf"
        File merged_bcf_index = "~{output_basename}.merged.bcf.csi"
    }
}

task BcftoolsConcatNaive {
    input {
        Array[File] bcfs
        String output_basename
    }

    Int disk_size_gb = ceil(2 * size(bcfs, "GiB")) + 5

    command <<<
        set -xeuo pipefail
        
        bcftools concat --naive -O b -o ~{output_basename}.bcf ~{sep=" " bcfs}
        bcftools index ~{output_basename}.bcf
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: "4 GiB"
        cpu: 1
        preemptible: 3
        noAddress: true
    }

    output {
        File output_bcf = "~{output_basename}.bcf"
        File output_bcf_index = "~{output_basename}.bcf.csi"
    }
}
