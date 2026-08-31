version 1.0

workflow Glimpse2LowPassImputationBatch {
    String pipeline_version = "1.1.0"

    input {
        Array[String] contigs
        Array[String] contig_groups # e.g., ["chr1,chr2,chr3", "chr4,chr5,chr6,chr7", ...] must cover all contigs continuously in order
        
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

    # Step 1: Pre-calculate reference files and counts
    scatter(contig_index in range(length(contigs))) {
        String current_contig = contigs[contig_index]
        Array[String] current_bcftools_shards = bcftools_shard_map[current_contig]
        Int current_shard_count = length(current_bcftools_shards)

        File sites_vcfs_tmp = reference_panel_prefix + "sites." + current_contig + ".vcf.gz"
        File sites_vcf_indices_tmp = reference_panel_prefix + "sites." + current_contig + ".vcf.gz.tbi"
        File sites_tables_tmp = reference_panel_prefix + "sites_table." + current_contig + ".gz"
        File sites_table_indices_tmp = reference_panel_prefix + "sites_table." + current_contig + ".gz.tbi"
        File reference_chunks_tmp = reference_panel_prefix + "reference_chunks." + current_contig + ".txt"

        call GetShards {
            input: 
                reference_chunks_file = reference_chunks_tmp
        }
    }

    # Step 2: Calculate global shard indices and group offsets in a single lightweight VM
    call CalculateShardIndices {
        input:
            shard_counts = current_shard_count
    }

    # Step 3: Extract utilizing explicit Contig Groups to balance preemption safety vs CRAM localization
    scatter(sample_idx in range(length(ParseCramManifest.crams))) {
        scatter(group_idx in range(length(contig_groups))) {
            
            call ExtractGenotypeLikelihoods {
                input:
                    cram = ParseCramManifest.crams[sample_idx],
                    cram_index = ParseCramManifest.cram_indices[sample_idx],
                    contig_group = contig_groups[group_idx],
                    contig_shard_starts = CalculateShardIndices.starts,
                    contigs = contigs,
                    flat_shards = flatten(current_bcftools_shards),
                    shard_counts = current_shard_count,
                    sites_vcfs = sites_vcfs_tmp,
                    sites_vcf_indices = sites_vcf_indices_tmp,
                    sites_tables = sites_tables_tmp,
                    sites_table_indices = sites_table_indices_tmp,
                    fasta = fasta,
                    fasta_index = fasta_index,
                    call_indels = call_indels
            }
        }
        
        # Flattens the scattered BCF outputs into a single sequential array for this sample
        Array[File] sample_flat_bcfs = flatten(ExtractGenotypeLikelihoods.flat_bcfs)
    }

    # WDL Native Transpose: Converts [Samples][Total_Shards] -> [Total_Shards][Samples] instantly
    Array[Array[File]] global_shard_x_sample = transpose(sample_flat_bcfs)

    # Step 4: Index Arithmetic loop to process the FLAT paste natively
    scatter(contig_index in range(length(contigs))) {
        String current_contig_for_merge = contigs[contig_index]
        Array[Int] global_indices_for_contig = CalculateShardIndices.indices[contig_index]

        scatter (local_shard_idx in range(length(global_indices_for_contig))) {
            Int global_shard_idx = global_indices_for_contig[local_shard_idx]
            Array[File] shard_sample_bcfs = global_shard_x_sample[global_shard_idx]
            
            # Execute flat paste directly
            call MergeSampleChunksBcfsWithPaste {
                input:
                    input_bcfs = shard_sample_bcfs,
                    output_basename = output_basename + "." + current_contig_for_merge + ".shard_" + local_shard_idx
            }
        }

        call BcftoolsConcatNaive {
            input:
                bcfs = MergeSampleChunksBcfsWithPaste.output_bcf,
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
        EOF
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk 10 HDD"
        memory: "1 GiB"
        preemptible: 3
        noAddress: true
    }

    output {
        Array[Array[Int]] indices = read_json("indices.json")
        Array[Int] starts = read_json("starts.json")
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
        disks: "local-disk 10 HDD"
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

task GetShards {
    input {
        File reference_chunks_file
    }

    command <<<
        python3 <<EOF
        import pandas as pd
        df = pd.read_csv('~{reference_chunks_file}', sep='\t', header=None, usecols=[0,1,2])
        df[1].to_csv('shards.txt', index=False, header=False)
        EOF
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        disks: "local-disk 10 HDD"
        memory: "1 GiB"
        preemptible: 3
        noAddress: true
    }

    output {
        Array[String] shards = read_lines("shards.txt")
    }
}

task ExtractGenotypeLikelihoods {
    input {
        File cram
        File cram_index
        
        String contig_group
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

        IFS=',' read -ra GROUP_CONTIGS <<< "~{contig_group}"
        
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
            
            vcf="vcf_${i}.vcf.gz"
            tbl="tbl_${i}.gz"
            count="${SHARD_COUNTS[$i]}"
            shard_idx="${SHARD_STARTS[$i]}"
            
            # Zero-padding ensures lexicographical glob matches the intended array sorting
            pad_contig=$(printf "%03d" $i)
            
            for (( j=0; j<count; j++ )); do
                shard="${ALL_SHARDS[$shard_idx]}"
                pad_shard=$(printf "%04d" $j)
                
                out_bcf="out_c${pad_contig}_s${pad_shard}.bcf"
                shard_tbl="tbl_c${pad_contig}_s${pad_shard}.tsv.gz"
                
                tabix -H "${tbl}" > "tmp.tsv" || true
                tabix "${tbl}" "${shard}" >> "tmp.tsv"
                bgzip -c "tmp.tsv" > "${shard_tbl}"
                tabix -s1 -b2 -e2 "${shard_tbl}"
                
                bcftools mpileup --no-version -f ~{fasta} ~{if !call_indels then "-I " else ""} --seed ~{seed} -E -a 'FORMAT/DP,FORMAT/AD' -r "${shard}" -T "${vcf}" -Ou ~{basename(cram)} | \
                bcftools call --no-version -Aim -C alleles -T "${shard_tbl}" -Ou | \
                bcftools norm --no-version -m -both -Ob -o "${out_bcf}"
                
                rm "tmp.tsv" "${shard_tbl}" "${shard_tbl}.tbi"
                
                ((++shard_idx))
            done
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
        # Lexicographical glob prevents PAPI indirect file array 404 bugs
        Array[File] flat_bcfs = glob("out_c*_s*.bcf")
    }
}

task MergeSampleChunksBcfsWithPaste {
    input {
        Array[File] input_bcfs
        String output_basename

        Int disk_size_gb = ceil(2.5 * size(input_bcfs, "GiB") + 20)
        Int mem_gb = 8
        Int cpu = 4
        Int preemptible = 1
    }

    command <<<
        set -euo pipefail

        bcfs=(~{sep=" " input_bcfs})

        mkfifo fifo_0
        mkfifo fifo_to_paste_0

        i=1

        fifos_to_paste=()
        md5sums=()
        
        bcftools view -h --no-version "${bcfs[0]}" | awk '!/^#CHROM/' > header.vcf

        (bcftools view -h --no-version "${bcfs[0]}" | grep '^#CHROM'; bcftools view -H "${bcfs[0]}") > fifo_0 &

        cat fifo_0 | tee fifo_to_paste_0 | cut -f1-5,9 | md5sum > md5sum_0 &

        for bcf in "${bcfs[@]:1}"; do
            fifo_name="fifo_$i"
            mkfifo "$fifo_name"

            fifo_name_to_md5="fifo_to_md5_$i"
            mkfifo "$fifo_name_to_md5"

            fifo_name_to_paste="fifo_to_paste_$i"
            mkfifo "$fifo_name_to_paste"
            fifos_to_paste+=("$fifo_name_to_paste")

            file_name_md5sum="md5sum_$i"
            md5sums+=("$file_name_md5sum")

            (bcftools view -h --no-version "${bcf}" | grep '^#CHROM'; bcftools view -H "${bcf}") > "$fifo_name" &
            cat "$fifo_name" | tee "$fifo_name_to_md5" | cut -f 10- > "$fifo_name_to_paste" &
            cut -f1-5,9 "$fifo_name_to_md5" | md5sum > "$file_name_md5sum" &

            ((i++))
        done

        mkfifo fifo_to_cat

        paste fifo_to_paste_0 ${fifos_to_paste[@]+"${fifos_to_paste[@]}"} | tee fifo_to_cat | awk 'NR % 5000000 == 0' | cut -f 1-5 &

        cat header.vcf fifo_to_cat | bcftools view --no-version -O b -o ~{output_basename}.bcf
        bcftools index ~{output_basename}.bcf

        for md5sum_file in ${md5sums[@]+"${md5sums[@]}"}; do
            diff <(cat md5sum_0) <(cat "$md5sum_file") >> /dev/null || (echo "Fields 1-5,9 do not match for $md5sum_file" && exit 1)
        done

        for fifo in fifo_*; do
            rm $fifo
        done
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
        disks: "local-disk " + disk_size_gb + " SSD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: 1
        noAddress: true
    }

    output {
        File output_bcf = "~{output_basename}.bcf"
        File output_bcf_index = "~{output_basename}.bcf.csi"
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
        disks: "local-disk " + disk_size_gb + " SSD"
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
