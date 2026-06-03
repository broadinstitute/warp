version 1.0

# This workflow performs low pass imputation using GLIMPSE2. It's designed to scale
# to approximately 1000 samples and be used as a subworkflow for Glimpse2LowPassImputation.wdl,
# which can handle larger sample sizes by splitting into batches and then merging results.

workflow Glimpse2LowPassImputationBatch {
    # if this changes, update the batch_pipeline_version value in Glimpse2LowPassImputation.wdl
    String pipeline_version = "0.0.8"

    input {

        Array[String] contigs

        # this is the path to a directory that contains sites vcf, sites table, and reference chunks file. should end with a "/"
        String reference_panel_prefix

        Array[File] crams
        Array[File] cram_indices
        Array[String] sample_ids
        File fasta
        File fasta_index
        String output_basename

        Boolean impute_reference_only_variants = false
        Boolean call_indels = false

        # batch size used when calling SplitIntoBatches to make variant calls from the crams
        Int calling_batch_size = 100

        # override for cpu used for glimpse phase task. Mostly used to set to 1 for determinism in testing
        Int? glimpse_phase_cpu_override

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.0.0"
        String glimpse_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2:1.0.0-2cee597-1778869818"
    }

    # we need to define this here so that it can be used in nested scatters below. Cromwell doesn't understand optional inputs
    # to tasks that are inside nested scatters, so we need to define a non-optional variable that we can use to pass the
    # value down to the GlimpsePhase task. If not defined, Cromwell fails the workflow
    Int defined_glimpse_phase_cpu_override = select_first([glimpse_phase_cpu_override, 4])

    if (length(crams) > 1) {
        call SplitIntoBatches {
            input:
                batch_size = calling_batch_size,
                crams = crams,
                cram_indices = cram_indices,
                sample_ids = sample_ids
        }
    }

    scatter(contig in contigs) {
        File sites_vcf = reference_panel_prefix + "sites." + contig + ".vcf.gz"
        File sites_vcf_index =reference_panel_prefix + "sites." + contig + ".vcf.gz.tbi"
        File sites_table = reference_panel_prefix + "sites_table." + contig + ".gz"
        File sites_table_index = reference_panel_prefix + "sites_table." + contig + ".gz.tbi"
        File reference_chunks = reference_panel_prefix + "reference_chunks." + contig + ".txt"

        Array[Array[String]] crams_batches = select_first([SplitIntoBatches.crams_batches, [select_first([crams])]])
        Array[Array[String]] cram_indices_batches = select_first([SplitIntoBatches.cram_indices_batches, [select_first([cram_indices])]])
        Array[Array[String]] sample_ids_batches = select_first([SplitIntoBatches.sample_ids_batches, [sample_ids]])

        scatter(i in range(length(crams_batches))) {
            call BcftoolsMpileup {
                input:
                    crams = crams_batches[i],
                    cram_indices = cram_indices_batches[i],
                    sample_ids = sample_ids_batches[i],
                    fasta = fasta,
                    fasta_index = fasta_index,
                    call_indels = call_indels,
                    sites_vcf = sites_vcf,
            }

            call BcftoolsCall {
                input:
                    mpileup_bcf = BcftoolsMpileup.output_bcf,
                    sites_table = sites_table,
                    sites_table_index = sites_table_index,
            }

            call BcftoolsNorm {
                input:
                    calls_bcf = BcftoolsCall.output_bcf,
            }
        }

        if (length(BcftoolsNorm.output_vcf) > 1) {
            call BcftoolsMerge {
                input:
                    vcfs = BcftoolsNorm.output_vcf,
                    vcf_indices = BcftoolsNorm.output_vcf_index,
                    output_basename = output_basename
            }
        }

        File phase_input_vcf = select_first([BcftoolsMerge.merged_vcf, BcftoolsNorm.output_vcf[0]])
        File phase_input_vcf_index = select_first([BcftoolsMerge.merged_vcf_index, BcftoolsNorm.output_vcf_index[0]])

        call ComputeShardsAndMemoryPerShard {
            input:
                reference_chunks_memory = reference_chunks
        }

        scatter (reference_chunk_index in range(length(ComputeShardsAndMemoryPerShard.reference_chunk_file_paths))) {

            call GlimpsePhase {
                input:
                    reference_chunk = ComputeShardsAndMemoryPerShard.reference_chunk_file_paths[reference_chunk_index],
                    input_vcf = phase_input_vcf,
                    input_vcf_index = phase_input_vcf_index,
                    impute_reference_only_variants = impute_reference_only_variants,
                    call_indels = call_indels,
                    sample_ids = sample_ids,
                    fasta = fasta,
                    fasta_index = fasta_index,
                    mem_gb = ComputeShardsAndMemoryPerShard.mem_gb_per_chunk[reference_chunk_index],
                    cpu = defined_glimpse_phase_cpu_override,
                    docker = glimpse_docker
            }
        }

        call GlimpseLigate {
            input:
                imputed_chunks = GlimpsePhase.imputed_vcf,
                imputed_chunks_indices = GlimpsePhase.imputed_vcf_index,
                output_basename = output_basename,
                docker = glimpse_docker
        }
    }

    output {
        Array[File] imputed_contig_ligated_vcfs = GlimpseLigate.imputed_vcf
        Array[File] imputed_contig_ligated_vcf_indices = GlimpseLigate.imputed_vcf_index
    }
}

task SplitIntoBatches {
    input {
        Int batch_size

        Array[String] crams
        Array[String] cram_indices
        Array[String] sample_ids
    }

    command <<<
        cat <<EOF > script.py
        import json

        batch_size = ~{batch_size}
        crams = ['~{sep="', '" crams}']
        cram_indices = ['~{sep="', '" cram_indices}']
        sample_ids = ['~{sep="', '" sample_ids}']

        crams_batches = [crams[i:i + batch_size] for i in range(0, len(crams), batch_size)]
        cram_indices_batches = [cram_indices[i:i + batch_size] for i in range(0, len(cram_indices), batch_size)]
        sample_ids_batches = [sample_ids[i:i + batch_size] for i in range(0, len(sample_ids), batch_size)]

        with open('crams.json', 'w') as json_file:
            json.dump(crams_batches, json_file)
        with open('cram_indices.json', 'w') as json_file:
            json.dump(cram_indices_batches, json_file)
        with open('sample_ids.json', 'w') as json_file:
            json.dump(sample_ids_batches, json_file)
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
        Array[Array[String]] crams_batches = read_json('crams.json')
        Array[Array[String]] cram_indices_batches = read_json('cram_indices.json')
        Array[Array[String]] sample_ids_batches = read_json('sample_ids.json')
    }
}

task ComputeShardsAndMemoryPerShard {
    input {
        File reference_chunks_memory
    }

    command <<<
        python3 << EOF
        import pandas as pd
        import numpy as np


        df = pd.read_csv('~{reference_chunks_memory}', sep='\t', header=None, usecols=[0,1,2], names=['contig', 'reference_shard', 'base_gb'])

        # write out reference shards to process
        df['reference_shard'].to_csv('reference_shard_file_paths.tsv', sep='\t', index=False, header=None)

        # calculate memory usage and save to file
        df['mem_gb'] = df['base_gb']
        df['mem_gb'] = df['mem_gb'].apply(lambda x: min(256, int(np.ceil(x))))  # cap at 256 GB
        df['mem_gb'].to_csv('memory_per_chunk.tsv', sep='\t', index=False, header=None)
        EOF
    >>>

    runtime {
        docker : "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        noAddress: true
    }

    output {
        Array[String] reference_chunk_file_paths = read_lines("reference_shard_file_paths.tsv")
        Array[Int] mem_gb_per_chunk = read_lines("memory_per_chunk.tsv")
    }
}

task BcftoolsMpileup {
    input {
        Array[File] crams
        Array[File] cram_indices
        File fasta
        File fasta_index
        Boolean call_indels
        Array[String] sample_ids

        File sites_vcf

        Int seed = 12345
        Int mem_gb = 6
        Int cpu = 1
        Int preemptible = 0
        Int max_retries = 3
    }

    Int disk_size_gb = ceil(1.5*size(crams, "GiB") + size(fasta, "GiB") + size(sites_vcf, "GiB")) + 10

    command <<<
        set -xeuo pipefail

        crams=(~{sep=' ' crams})
        sample_ids=(~{sep=' ' sample_ids})

        for i in "${!crams[@]}"; do
            echo "* ${crams[$i]} ${sample_ids[$i]}" >> sample_name_mapping.txt
        done

        bcftools mpileup -f ~{fasta} ~{if !call_indels then "-I" else ""} -G sample_name_mapping.txt --seed ~{seed} -E -a 'FORMAT/DP,FORMAT/AD' -T ~{sites_vcf} -Ou -o mpileup.bcf ~{sep=" " crams}
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
        noAddress: true
    }

    output {
        File output_bcf = "mpileup.bcf"
    }
}

task BcftoolsCall {
    input {
        File mpileup_bcf

        File sites_table
        File sites_table_index

        Int mem_gb = 12
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 3
    }

    Int disk_size_gb = ceil(3*size(mpileup_bcf, "GiB") + size(sites_table, "GiB")) + 10

    command <<<
        set -xeuo pipefail

        bcftools call -Aim -C alleles -T ~{sites_table} -Ou ~{mpileup_bcf} -o calls.bcf
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
        disks: "local-disk " + disk_size_gb + " SSD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
        noAddress: true
    }

    output {
        File output_bcf = "calls.bcf"
    }
}

task BcftoolsNorm {
    input {
        File calls_bcf

        Int mem_gb = 6
        Int cpu = 1
        Int preemptible = 3
        Int max_retries = 3
    }

    Int disk_size_gb = ceil(3*size(calls_bcf, "GiB")) + 10

    command <<<
        set -xeuo pipefail


        bcftools norm -m -both -Oz -o normalized.vcf.gz ~{calls_bcf}
        bcftools index -t normalized.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
        disks: "local-disk " + disk_size_gb + " SSD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
        noAddress: true
    }

    output {
        File output_vcf = "normalized.vcf.gz"
        File output_vcf_index = "normalized.vcf.gz.tbi"
    }
}

task BcftoolsMerge {
    input {
        Array[File] vcfs
        Array[File] vcf_indices
        Int mem_gb = 6
        Int cpu = 1
        Int preemptible = 0
        Int max_retries = 3

        String output_basename
    }

    Int disk_size_gb = ceil(3*size(vcfs, "GiB")) + 50

    command <<<
        set -euo pipefail
        bcftools merge -O z -o ~{output_basename}.bcftools.merged.vcf.gz ~{sep=" " vcfs}
        bcftools index -t ~{output_basename}.bcftools.merged.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
        noAddress: true
    }

    output {
        File merged_vcf = "~{output_basename}.bcftools.merged.vcf.gz"
        File merged_vcf_index = "~{output_basename}.bcftools.merged.vcf.gz.tbi"
    }
}

task GlimpsePhase {
    input {
        File? input_vcf
        File? input_vcf_index
        Array[File]? crams
        Array[File]? cram_indices
        Array[String] sample_ids
        File? fasta
        File? fasta_index
        File reference_chunk

        Boolean impute_reference_only_variants
        Boolean call_indels
        Int? n_burnin
        Int? n_main
        Int? effective_population_size
        Int seed = 15052011

        Int cpu = 4 # note that setting cpu > 1 will introduce non-determinism in GLIMPSE Phase due to multi-threading
        Int mem_gb = 16
        Int disk_size_gb = ceil(2.2 * size(input_vcf, "GiB") + size(reference_chunk, "GiB") + 0.003 * length(select_first([crams, []])) + 10)
        Int preemptible = 30
        Int max_retries = 3
        String docker
    }

    parameter_meta {
        crams: {
                   localization_optional: true
               }
        cram_indices: {
                          localization_optional: true
                      }
        input_vcf: {
                   localization_optional: true
               }
        input_vcf_index: {
                          localization_optional: true
                      }
    }

    String bam_file_list_input = if defined(crams) then "--bam-list crams.list" else ""
    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(/google-cloud-sdk/bin/gcloud auth application-default print-access-token)

        cram_paths=( ~{sep=" " crams} )
        cram_index_paths=( ~{sep=" " cram_indices} )
        sample_ids=( ~{sep=" " sample_ids} )

        duplicate_cram_filenames=$(printf "%s\n" "${cram_paths[@]}" | xargs -I {} basename {} | sort | uniq -d)
        if [ ! -z "$duplicate_cram_filenames" ]; then
            echo "ERROR: The input CRAMs contain multiple files with the same basename, which leads to an error due to the way that htslib is implemented. Duplicate filenames:"
            printf "%s\n" "${duplicate_cram_filenames[@]}"
            exit 1
        fi

        if ~{if defined(cram_indices) then "true" else "false"}; then
            for i in "${!cram_paths[@]}" ; do
                echo -e "${cram_paths[$i]}##idx##${cram_index_paths[$i]} ${sample_ids[$i]}" >> crams.list
            done
        else
            for i in "${!cram_paths[@]}"; do
                echo -e "${cram_paths[$i]} ${sample_ids[$i]}" >> crams.list
            done
        fi

        cmd="/bin/GLIMPSE2_phase \
        ~{"--input-gl " + input_vcf} \
        --reference ~{reference_chunk} \
        --output phase_output.bcf \
        --threads ~{cpu} \
        --seed ~{seed} \
        ~{if impute_reference_only_variants then "--impute-reference-only-variants" else ""} ~{if call_indels then "--call-indels" else ""} \
        ~{"--burnin " + n_burnin} ~{"--main " + n_main} \
        ~{"--ne " + effective_population_size} \
        ~{bam_file_list_input} \
        ~{"--fasta " + fasta} \
        --checkpoint-file-out checkpoint.bin"

        if [ -s "checkpoint.bin" ]; then
            cmd="$cmd --checkpoint-file-in checkpoint.bin"
        fi



        #check for read error which corresponds exactly to end of cram/bam block.
        #This currently triggers a warning message from htslib, but doesn't return any error.
        #We need to make sure that stderr is maintained since cromwell looks for oom strings
        #in stderr

        eval $cmd 2> >(tee glimpse_stderr.log >&2)

        if grep -q "EOF marker is absent" glimpse_stderr.log; then
            echo "An input file appears to be truncated.  This may be either a truly truncated file which needs to be fixed, or a networking error which can just be retried."
        exit 1
        fi
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " SSD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
        checkpointFile: "checkpoint.bin"
        noAddress: true
    }

    output {
        File imputed_vcf = "phase_output.bcf"
        File imputed_vcf_index = "phase_output.bcf.csi"
    }
}

task GlimpseLigate {
    input {
        Array[File] imputed_chunks
        Array[File] imputed_chunks_indices
        String output_basename

        Int mem_gb = 4
        Int cpu = 2
        Int disk_size_gb = ceil(3 * size(imputed_chunks, "GiB") + 100)
        Int preemptible = 0
        Int max_retries = 3
        String docker
    }

    command <<<
        set -xeuo pipefail

        NPROC=$(nproc)
        echo "nproc reported ${NPROC} CPUs, using that number as the threads argument for GLIMPSE."

        /bin/GLIMPSE2_ligate --input ~{write_lines(imputed_chunks)} --output ~{output_basename}.imputed.vcf.gz --threads ${NPROC}

        # GLIMPSE2_ligate creates an index, but it is not compatible with GATK tools so we regenerate it with tabix
        tabix -f ~{output_basename}.imputed.vcf.gz
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
        noAddress: true
    }

    output {
        File imputed_vcf = "~{output_basename}.imputed.vcf.gz"
        File imputed_vcf_index = "~{output_basename}.imputed.vcf.gz.tbi"
    }
}
