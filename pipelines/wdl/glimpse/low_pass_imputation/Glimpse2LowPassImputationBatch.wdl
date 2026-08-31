version 1.0

import "../../../../tasks/wdl/Glimpse2LowPassImputationTasks.wdl" as Glimpse2LowPassImputationTasks

workflow Glimpse2LowPassImputationBatch {
    String pipeline_version = "1.1.0"

    input {
        Array[String] contigs
        String reference_panel_prefix

        File cram_manifest
        File fasta
        File fasta_index
        String output_basename

        Boolean impute_reference_only_variants = false
        Boolean call_indels = false

        Map[String, Array[String]] bcftools_shard_map
        Array[Int] hierarchical_merge_batch_sizes

        Int? glimpse_phase_cpu_override

        String gatk_docker
        String glimpse_docker
    }

    Int defined_glimpse_phase_cpu_override = select_first([glimpse_phase_cpu_override, 4])

    call ParseCramManifest {
        input:
            cram_manifest = cram_manifest
    }

    # Step 1: Pre-calculate reference files and flattened shard counts
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

    # Step 2: CRAMs are localized exactly ONCE per sample. 
    scatter(sample_idx in range(length(ParseCramManifest.crams))) {
        call ExtractGenotypeLikelihoods {
            input:
                cram = ParseCramManifest.crams[sample_idx],
                cram_index = ParseCramManifest.cram_indices[sample_idx],
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

    # Step 3: Loop through contigs to perform the hierarchical merge
    scatter(contig_index in range(length(contigs))) {
        
        String current_contig_for_merge = contigs[contig_index]
        Array[String] current_bcftools_shards_for_merge = bcftools_shard_map[current_contig_for_merge]

        scatter(sample_idx in range(length(ParseCramManifest.crams))) {
            Array[File] sample_shards_for_contig = ExtractGenotypeLikelihoods.out_bcfs[sample_idx][contig_index]
        }
        
        Array[Array[File]] shard_sample_bcfs = transpose(sample_shards_for_contig)

        scatter (shard_idx in range(length(current_bcftools_shards_for_merge))) {
            
            # Level 1
            call ChunkBcfs as Level1Chunk {
                input:
                    bcfs = shard_sample_bcfs[shard_idx],
                    batch_size = hierarchical_merge_batch_sizes[0]
            }
            scatter (l1_idx in range(length(Level1Chunk.out_chunks))) {
                call MergeSampleChunksBcfsWithPaste as Level1Merge {
                    input:
                        input_bcfs = Level1Chunk.out_chunks[l1_idx],
                        output_basename = output_basename + "." + current_contig_for_merge + ".shard_" + shard_idx + ".l1_" + l1_idx
                }
            }

            # Level 2 (Conditional)
            if (length(Level1Merge.output_bcf) > 1) {
                call ChunkBcfs as Level2Chunk {
                    input:
                        bcfs = Level1Merge.output_bcf,
                        batch_size = hierarchical_merge_batch_sizes[1]
                }
                scatter (l2_idx in range(length(Level2Chunk.out_chunks))) {
                    call MergeSampleChunksBcfsWithPaste as Level2Merge {
                        input:
                            input_bcfs = Level2Chunk.out_chunks[l2_idx],
                            output_basename = output_basename + "." + current_contig_for_merge + ".shard_" + shard_idx + ".l2_" + l2_idx
                    }
                }
            }
            
            # Select the final shard depending on whether Level 2 ran
            Array[File] resolved_final_shard_array = select_first([Level2Merge.output_bcf, Level1Merge.output_bcf])
            File final_shard_bcf = resolved_final_shard_array[0]
        }

        call BcftoolsConcat {
            input:
                bcfs = final_shard_bcf,
                output_basename = output_basename + "." + current_contig_for_merge + ".concat"
        }
    }

    output {
        Int total_samples = ParseCramManifest.total_samples
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

        # Read the manifest
        df = pd.read_csv("~{cram_manifest}", sep='\t')

        # Check for required columns
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

        # Co-locate the CRAM and index in the local task directory
        ln -sf ~{cram} "~{basename(cram)}"
        ln -sf ~{cram_index} "~{basename(cram)}.crai"

        # Load WDL arrays into Bash parallel arrays
        CONTIGS=(~{sep=" " contigs})
        VCFS=(~{sep=" " sites_vcfs})
        TABLES=(~{sep=" " sites_tables})
        SHARD_COUNTS=(~{sep=" " shard_counts})

        # Load the flattened shards array line-by-line safely
        ALL_SHARDS=()
        while IFS= read -r line; do
            ALL_SHARDS+=("$line")
        done < "~{write_lines(flat_shards)}"

        echo "[" > outputs.json

        shard_idx=0
        for i in "${!CONTIGS[@]}"; do
            chrom="${CONTIGS[$i]}"
            vcf="${VCFS[$i]}"
            tbl="${TABLES[$i]}"
            count="${SHARD_COUNTS[$i]}"
            
            if [ $i -gt 0 ]; then
                echo "," >> outputs.json
            fi
            echo "  [" >> outputs.json
            
            for (( j=0; j<count; j++ )); do
                shard="${ALL_SHARDS[$shard_idx]}"
                
                # Format string for safe filenames
                safe_shard=$(echo "$shard" | tr ':' '_' | tr '-' '_')
                out_bcf="output.${chrom}.${safe_shard}.bcf"
                
                # Execute bcftools piped pipeline
                bcftools mpileup -f ~{fasta} ~{if !call_indels then "-I " else ""} --seed ~{seed} -E -a 'FORMAT/DP,FORMAT/AD' -r "${shard}" -T "${vcf}" -Ou ~{basename(cram)} | \
                bcftools call -Aim -C alleles -T "${tbl}" -Ou | \
                bcftools norm -m -both -Ob -o "${out_bcf}"
                
                bcftools index "${out_bcf}"
                
                # Write to the 2D JSON array file for WDL to read back
                if [ $j -gt 0 ]; then
                    echo "    ,\"${out_bcf}\"" >> outputs.json
                else
                    echo "    \"${out_bcf}\"" >> outputs.json
                fi
                
                ((shard_idx++))
            done
            
            echo "  ]" >> outputs.json
        done

        echo "]" >> outputs.json
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
        Array[Array[File]] out_bcfs = read_json("outputs.json")
    }
}

task ChunkBcfs {
    input {
        Array[File] bcfs
        Int batch_size
    }

    command <<<
        python3 <<EOF
        import json
        bcfs = ['~{sep="', '" bcfs}']
        bs = ~{batch_size}
        chunks = [bcfs[i:i + bs] for i in range(0, len(bcfs), bs)]
        with open("chunks.json", "w") as f:
            json.dump(chunks, f)
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
        Array[Array[File]] out_chunks = read_json("chunks.json")
    }
}

task MergeSampleChunksBcfsWithPaste {
    input {
        Array[File] input_bcfs
        String output_basename

        Int disk_size_gb = ceil(2.5 * size(input_bcfs, "GiB") + 20)
        Int mem_gb = 8
        Int cpu = 4
        Int preemptible = 0
    }

    command <<<
        set -euo pipefail

        bcfs=(~{sep=" " input_bcfs})

        mkfifo fifo_0
        mkfifo fifo_to_paste_0

        i=1

        fifos_to_paste=()
        md5sums=()
        
        # 1. Extract header safely without the #CHROM line
        bcftools view -h ${bcfs[0]} | awk '!/^#CHROM/' > header.vcf

        # 2. Extract #CHROM line and raw data lines for the first BCF
        (bcftools view -h ${bcfs[0]} | grep '^#CHROM'; bcftools view -H ${bcfs[0]}) > fifo_0 &

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

            # Extract #CHROM line and raw data lines for sequential BCFs
            (bcftools view -h ${bcf} | grep '^#CHROM'; bcftools view -H ${bcf}) > "$fifo_name" &
            cat "$fifo_name" | tee "$fifo_name_to_md5" | cut -f 10- > "$fifo_name_to_paste" &
            cut -f1-5,9 "$fifo_name_to_md5" | md5sum > "$file_name_md5sum" &

            ((i++))
        done

        mkfifo fifo_to_cat

        paste fifo_to_paste_0 "${fifos_to_paste[@]}" | tee fifo_to_cat | awk 'NR % 5000000 == 0' | cut -f 1-5 &

        # 3. Stitch header and pasted text back together, and output directly to compressed BCF
        cat header.vcf fifo_to_cat | bcftools view -O b -o ~{output_basename}.bcf
        bcftools index ~{output_basename}.bcf

        for md5sum_file in "${md5sums[@]}"; do
            diff <(cat md5sum_0) <(cat $md5sum_file) >> /dev/null || (echo "Fields 1-5,9 do not match for $md5sum_file" && exit 1)
        done

        for fifo in fifo_*; do
            rm $fifo
        done
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
        disks: "local-disk " + disk_size_gb + " HDD"
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

task BcftoolsConcat {
    input {
        Array[File] bcfs
        String output_basename
    }

    Int disk_size_gb = ceil(2 * size(bcfs, "GiB") + 5)

    command <<<
        set -xeuo pipefail
        
        bcftools concat -O b -o ~{output_basename}.bcf ~{sep=" " bcfs} --naive
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
