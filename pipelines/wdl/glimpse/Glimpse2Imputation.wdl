version 1.0

workflow Glimpse2Imputation {
    input {
        # List of files, one per line
        File reference_chunks
        File sites_vcf

        Int bcftools_threads
        Int calling_batch_size
        Int calling_mem_gb

        File? input_vcf
        File? input_vcf_index
        Array[File]? crams
        Array[File]? cram_indices
        Array[String] sample_ids
        File fasta
        File fasta_index
        String output_basename

        File ref_dict

        Boolean impute_reference_only_variants = false
        Boolean call_indels = false
        Int? n_burnin
        Int? n_main
        Int? effective_population_size

        Int preemptible = 9
        String docker = "us.gcr.io/broad-dsde-methods/glimpse:kachulis_ck_bam_reader_retry_cf5822c"
        String docker_extract_num_sites_from_reference_chunk = "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a"
        Int cpu_ligate = 4
        Int mem_gb_ligate = 4
        Int? cpu_phase
        Int? mem_gb_phase
        File? monitoring_script
    }

    if (defined(crams)) {
        if (length(select_first([crams])) > 1) {
            call SplitIntoBatches {
                input:
                    batch_size = calling_batch_size,
                    crams = select_first([crams]),
                    cram_indices = select_first([cram_indices]),
                    sample_ids = sample_ids
            }
        }
        Array[Array[String]] crams_batches = select_first([SplitIntoBatches.crams_batches, [select_first([crams])]])
        Array[Array[String]] cram_indices_batches = select_first([SplitIntoBatches.cram_indices_batches, [select_first([cram_indices])]])
        Array[Array[String]] sample_ids_batches = select_first([SplitIntoBatches.sample_ids_batches, [select_first([sample_ids])]])

        scatter(i in range(length(crams_batches))) {
            call BcftoolsCall {
                input:
                    crams = crams_batches[i],
                    cram_indices = cram_indices_batches[i],
                    sample_ids = sample_ids_batches[i],
                    fasta = fasta,
                    fasta_index = fasta_index,
                    call_indels = call_indels,
                    sites_vcf = sites_vcf,
                    cpu = bcftools_threads,
                    mem_gb = calling_mem_gb
            }
        }

        if (length(BcftoolsCall.output_vcf) > 1) {
            call BcftoolsMerge {
                input:
                    vcfs = BcftoolsCall.output_vcf,
                    vcf_indices = BcftoolsCall.output_vcf_index,
                    output_basename = output_basename
            }
        }

        File merged_vcf = select_first([BcftoolsMerge.merged_vcf, BcftoolsCall.output_vcf[0]])
        File merged_vcf_index = select_first([BcftoolsMerge.merged_vcf_index, BcftoolsCall.output_vcf_index[0]])
    }

    scatter (reference_chunk in read_lines(reference_chunks)) {
        if (!defined(cpu_phase) || !defined(mem_gb_phase)) {
            call GetNumberOfSitesInChunk {
                input:
                    reference_chunk = reference_chunk,
                    docker = docker_extract_num_sites_from_reference_chunk
            }

            Int n_rare = GetNumberOfSitesInChunk.n_rare
            Int n_common = GetNumberOfSitesInChunk.n_common

            if (defined(input_vcf)) {
                call CountSamples {
                    input:
                        vcf = select_first([input_vcf])
                }
            }

            Int n_samples = select_first([CountSamples.nSamples, length(select_first([crams]))])

            call SelectResourceParameters {
                input:
                    n_rare = n_rare,
                    n_common = n_common,
                    n_samples = n_samples
            }

            if (SelectResourceParameters.memory_gb > 256 || SelectResourceParameters.request_n_cpus > 32) {
                # force failure if we're accidently going to request too much resources and spend too much money
                Int safety_check_memory_gb = -1
                Int safety_check_n_cpu = -1
            }
        }

        call GlimpsePhase {
            input:
                reference_chunk = reference_chunk,
                input_vcf = select_first([merged_vcf,input_vcf]),
                input_vcf_index = select_first([merged_vcf_index,input_vcf_index]),
                impute_reference_only_variants = impute_reference_only_variants,
                n_burnin = n_burnin,
                n_main = n_main,
                effective_population_size = effective_population_size,
                call_indels = call_indels,
                sample_ids = sample_ids,
                fasta = fasta,
                fasta_index = fasta_index,
                preemptible = preemptible,
                docker = docker,
                cpu = select_first([cpu_phase, safety_check_n_cpu, SelectResourceParameters.request_n_cpus]),
                mem_gb = select_first([mem_gb_phase, safety_check_memory_gb, SelectResourceParameters.memory_gb]),
                monitoring_script = monitoring_script
        }
    }

    call GlimpseLigate {
        input:
            imputed_chunks = GlimpsePhase.imputed_vcf,
            imputed_chunks_indices = GlimpsePhase.imputed_vcf_index,
            output_basename = output_basename,
            ref_dict = ref_dict,
            preemptible = preemptible,
            docker = docker,
            cpu = cpu_ligate,
            mem_gb = mem_gb_ligate,
    }

    if (length(select_all(GlimpsePhase.coverage_metrics)) > 0) {
        call CombineCoverageMetrics {
            input:
                cov_metrics = select_all(GlimpsePhase.coverage_metrics),
                output_basename = output_basename
        }
    }

    call CollectQCMetrics {
        input:
            imputed_vcf = GlimpseLigate.imputed_vcf,
            output_basename = output_basename,
            monitoring_script = monitoring_script
    }


    output {
        File imputed_vcf = GlimpseLigate.imputed_vcf
        File imputed_vcf_index = GlimpseLigate.imputed_vcf_index
        File imputed_vcf_md5sum = GlimpseLigate.imputed_vcf_md5sum

        File qc_metrics = CollectQCMetrics.qc_metrics
        File? coverage_metrics = CombineCoverageMetrics.coverage_metrics

        Array[File?] glimpse_phase_monitoring = GlimpsePhase.monitoring
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
    }

    output {
        Array[Array[String]] crams_batches = read_json('crams.json')
        Array[Array[String]] cram_indices_batches = read_json('cram_indices.json')
        Array[Array[String]] sample_ids_batches = read_json('sample_ids.json')
    }
}

task BcftoolsCall {
    input {
        Array[File] crams
        Array[File] cram_indices
        File fasta
        File fasta_index
        Boolean call_indels
        Array[String] sample_ids

        File sites_vcf
        File sites_table
        File sites_table_index

        Int mem_gb = 4
        Int cpu = 2
        Int preemptible = 3
    }

    Int disk_size_gb = ceil(1.5*size(crams, "GiB") + size(fasta, "GiB") + size(sites_table, "GiB")) + 10

    String out_basename = "batch"

    command <<<
        set -xeuo pipefail

        crams=(~{sep=' ' crams})
        sample_ids=(~{sep=' ' sample_ids})

        for i in "${!crams[@]}"; do
        echo "* ${crams[$i]} ${sample_ids[$i]}" >> sample_name_mapping.txt
        done

        bcftools mpileup -f ~{fasta} ~{if !call_indels then "-I" else ""} -G sample_name_mapping.txt -E -a 'FORMAT/DP,FORMAT/AD' -T ~{sites_vcf} -Ou ~{sep=" " crams} \
        | bcftools call -Aim -C alleles -T ~{sites_table} -Ou - \
        | bcftools norm -m -both -Oz -o ~{out_basename}.vcf.gz -
        bcftools index -t ~{out_basename}.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File output_vcf = "~{out_basename}.vcf.gz"
        File output_vcf_index = "~{out_basename}.vcf.gz.tbi"
    }
}

task BcftoolsMerge {
    input {
        Array[File] vcfs
        Array[File] vcf_indices
        Int mem_gb = 4
        Int cpu = 2
        Int preemptible = 1

        String output_basename
    }

    Int disk_size_gb = ceil(2*size(vcfs, "GiB")) + 10

    command <<<
        set -euo pipefail
        bcftools merge -O z -o ~{output_basename}.bcftools.merged.vcf.gz ~{sep=" " vcfs}
        bcftools index -t ~{output_basename}.bcftools.merged.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools:v1.3"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
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

        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(2.2 * size(input_vcf, "GiB") + size(reference_chunk, "GiB") + 0.003 * length(select_first([crams, []])) + 10)
        Int preemptible = 9
        Int max_retries = 3
        String docker
        File? monitoring_script
    }

    parameter_meta {
        crams: {
                   localization_optional: true
               }
        cram_indices: {
                          localization_optional: true
                      }
    }

    String bam_file_list_input = if defined(crams) then "--bam-list crams.list" else ""
    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(/root/google-cloud-sdk/bin/gcloud auth application-default print-access-token)
        ~{"bash " + monitoring_script + " > monitoring.log &"}

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
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
        checkpointFile: "checkpoint.bin"
    }

    output {
        File imputed_vcf = "phase_output.bcf"
        File imputed_vcf_index = "phase_output.bcf.csi"
        File? monitoring = "monitoring.log"
        File? coverage_metrics = "phase_output_stats_coverage.txt.gz"
    }
}

task GlimpseLigate {
    input {
        Array[File] imputed_chunks
        Array[File] imputed_chunks_indices
        String output_basename
        File ref_dict

        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(2.2 * size(imputed_chunks, "GiB") + 100)
        Int preemptible = 1
        Int max_retries = 3
        String docker
    }

    command <<<
        set -xeuo pipefail

        NPROC=$(nproc)
        echo "nproc reported ${NPROC} CPUs, using that number as the threads argument for GLIMPSE."

        /bin/GLIMPSE2_ligate --input ~{write_lines(imputed_chunks)} --output ligated.vcf.gz --threads ${NPROC}

        # Set correct reference dictionary
        bcftools view -h --no-version ligated.vcf.gz > old_header.vcf
        java -jar /picard.jar UpdateVcfSequenceDictionary -I old_header.vcf --SD ~{ref_dict} -O new_header.vcf
        bcftools reheader -h new_header.vcf -o ~{output_basename}.imputed.vcf.gz ligated.vcf.gz
        tabix ~{output_basename}.imputed.vcf.gz

        md5sum ~{output_basename}.imputed.vcf.gz | awk '{ print $1 }' > ~{output_basename}.imputed.vcf.gz.md5sum
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
    }

    output {
        File imputed_vcf = "~{output_basename}.imputed.vcf.gz"
        File imputed_vcf_index = "~{output_basename}.imputed.vcf.gz.tbi"
        File imputed_vcf_md5sum = "~{output_basename}.imputed.vcf.gz.md5sum"
    }
}

task CollectQCMetrics {
    input {
        File imputed_vcf
        String output_basename

        Int preemptible = 1
        String docker = "hailgenetics/hail:0.2.126-py3.11"
        Int cpu = 4
        Int mem_gb = 16
        File? monitoring_script
    }

    parameter_meta {
        imputed_vcf: {
                         localization_optional: true
                     }
    }

    Int disk_size_gb = 100

    command <<<
        set -euo pipefail

        ~{"bash " + monitoring_script + " > monitoring.log &"}

        cat <<'EOF' > script.py
        import hail as hl
        import pandas as pd

        # Calculate metrics
        hl.init(default_reference='GRCh38', idempotent=True)
        vcf = hl.import_vcf('~{imputed_vcf}', force_bgz=True)
        qc = hl.sample_qc(vcf)
        qc_pd = qc.cols().flatten() \
        .rename({'sample_qc.' + col: col for col in list(qc['sample_qc'])}) \
        .rename({'s': 'sample_id'}) \
        .to_pandas()
        qc_pd.to_csv('~{output_basename}.qc_metrics.tsv', sep='\t', index=False, float_format='%.4f')
        EOF
        python3 script.py
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File qc_metrics = "~{output_basename}.qc_metrics.tsv"
        File? monitoring = "monitoring.log"
    }
}

task GetNumberOfSitesInChunk {
    input {
        File reference_chunk

        String docker
        Int mem_gb = 4
        Int cpu = 4
        Int disk_size_gb = ceil(size(reference_chunk, "GiB") + 10)
        Int preemptible = 1
        Int max_retries = 3
    }

    command <<<
        set -xeuo pipefail
        /bin/GLIMPSE2_extract_num_sites_from_reference_chunk ~{reference_chunk} > n_sites.txt
        cat n_sites.txt
        grep "Lrare" n_sites.txt | sed 's/Lrare=//' > n_rare.txt
        grep "Lcommon" n_sites.txt | sed 's/Lcommon=//' > n_common.txt
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: max_retries
    }

    output {
        Int n_rare = read_int("n_rare.txt")
        Int n_common = read_int("n_common.txt")
    }
}

task CountSamples {
    input {
        File vcf

        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int cpu = 1
        Int memory_mb = 3000
        Int disk_size_gb = 10 + ceil(size(vcf, "GiB"))
    }

    command <<<
        bcftools query -l ~{vcf} | wc -l
    >>>

    runtime {
        docker: bcftools_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }
    output {
        Int nSamples = read_int(stdout())
    }
}

task SelectResourceParameters {
    input {
        Int n_rare
        Int n_common
        Int n_samples
    }

    command <<<
        python3 << EOF
        import math
        n_rare = ~{n_rare}
        n_common = ~{n_common}
        n_samples = ~{n_samples}
        n_sites = n_common + n_rare

        # try to keep expected runtime under 4 hours, but don't ask for more than 32 cpus, or 256 GB memory
        estimated_needed_threads = min(math.ceil(5e-6*n_sites*n_samples/240), 32)
        estimated_needed_memory_gb = min(math.ceil((800e-3 + 0.97e-6 * n_rare * estimated_needed_threads + 14.6e-6 * n_common * estimated_needed_threads + 6.5e-9 * (n_rare + n_common) * n_samples + 13.7e-3 * n_samples + 1.8e-6*(n_rare + n_common)*math.log(n_samples))), 256)
        # recalc allowable threads, may be some additional threads available due to rounding memory up
        threads_to_use = max(math.floor((estimated_needed_memory_gb - (800e-3 + 6.5e-9 * (n_rare + n_common) * n_samples + 13.7e-3 * n_samples + 1.8e-6*(n_rare + n_common)*math.log(n_samples)))/(0.97e-6 * n_rare + 14.6e-6 * n_common)), 1)
        #estimated_needed_memory_gb = math.ceil(1.2 * estimated_needed_memory_gb)

        with open("n_cpus_request.txt", "w") as f_cpus_request:
        f_cpus_request.write(f'{int(threads_to_use)}')

        with open("memory_gb.txt", "w") as f_mem:
        f_mem.write(f'{int(estimated_needed_memory_gb)}')
        EOF
    >>>

    runtime {
        docker : "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
    }

    output {
        Int memory_gb = read_int("memory_gb.txt")
        Int request_n_cpus = read_int("n_cpus_request.txt")
    }
}

task CombineCoverageMetrics
{
    input {
        Array[File] cov_metrics
        String output_basename
    }

    command <<<
        set -euo pipefail

        cov_files=( ~{sep=" " cov_metrics} )

        for i in "${!cov_files[@]}"; do
        if [ $i -eq 0 ]; then
        n_skip=1
        echo 'Chunk' > chunk_col.txt
        else
        n_skip=2
        fi
        # glimpse coverage metrics are formatted to be human readable in a command line, not machine readable or consistent.  ie, number of tabs
        # are variable between columns depending on length of sample names, odd things like that.  We want these to be machine readable tables,
        # so need to fix this.
        zcat ${cov_files[$i]} | tail -n +$((n_skip + 1)) | sed s/%//g | sed s/"No data"/"No data pct"/g | sed s/\\t\\t/\\t/g >> cov_file.txt
        n_lines_cov=$(< cov_file.txt wc -l)
        n_lines_chunk=$(< chunk_col.txt wc -l)
        n_lines_out=$((n_lines_cov-n_lines_chunk))
        echo 'n_lines_out=' ${n_lines_out}
        echo ${cov_files[$i]}
        { yes ${i} || :; } | head -n ${n_lines_out} >> chunk_col.txt
        done

        paste chunk_col.txt cov_file.txt > ~{output_basename}.coverage_metrics.txt

    >>>

    runtime {
        docker: "ubuntu:24.04"
    }

    output {
        File coverage_metrics="~{output_basename}.coverage_metrics.txt"
    }
}
