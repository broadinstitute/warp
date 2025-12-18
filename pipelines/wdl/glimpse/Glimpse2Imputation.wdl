version 1.0

workflow Glimpse2Imputation {
    input {
        # List of files, one per line
        File reference_chunks_memory    # File with baseline GB and slope per sample for memory calculation, for each chunk

        File? input_vcf
        File? input_vcf_index
        Array[File]? crams
        Array[File]? cram_indices
        Array[String] sample_ids
        File? fasta
        File? fasta_index

        String output_basename
        Array[String] contigs
        File ref_dict


        Boolean impute_reference_only_variants = false
        Boolean call_indels = false
        Int? n_burnin
        Int? n_main
        Int? effective_population_size
        
        Int preemptible = 9
        String docker = "us.gcr.io/broad-dsde-methods/glimpse:odelaneau_bd93ade"
        String docker_extract_num_sites_from_reference_chunk = "us.gcr.io/broad-dsde-methods/glimpse_extract_num_sites_from_reference_chunks:michaelgatzen_edc7f3a"
        Int cpu_ligate = 4
        Int mem_gb_ligate = 4
        Int? cpu_phase
        Int? mem_gb_phase
        Float extra_mem_scaling_phase = 0.0
    }

    if (defined(input_vcf)) {
        call CountSamples {
            input:
                vcf = select_first([input_vcf])
        }
    }

    Int n_samples = select_first([CountSamples.nSamples, length(select_first([crams]))])

    call ComputeShardsAndMemoryPerShard {
        input:
            reference_chunks_memory = reference_chunks_memory,
            contigs = contigs,
            n_samples = n_samples
    }

    scatter (reference_chunk in zip(ComputeShardsAndMemoryPerShard.reference_chunk_file_paths, ComputeShardsAndMemoryPerShard.mem_gb_per_chunk)) {
        call GlimpsePhase {
            input:
                reference_chunk = reference_chunk.left,
                input_vcf = input_vcf,
                input_vcf_index = input_vcf_index,
                impute_reference_only_variants = impute_reference_only_variants,
                n_burnin = n_burnin,
                n_main = n_main,
                effective_population_size = effective_population_size,
                call_indels = call_indels,
                crams = crams,
                cram_indices = cram_indices,
                sample_ids = sample_ids,
                fasta = fasta,
                fasta_index = fasta_index,
                preemptible = preemptible,
                docker = docker,
                cpu = select_first([cpu_phase, 1]),
                mem_gb = select_first([mem_gb_phase, reference_chunk.right])
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

    call CombineCoverageMetrics {
        input:
            cov_metrics = GlimpsePhase.coverage_metrics,
            output_basename = output_basename
    }
   
    call CollectQCMetrics {
        input:
            imputed_vcf = GlimpseLigate.imputed_vcf,
            output_basename = output_basename,
    }
    

    output {
        File imputed_vcf = GlimpseLigate.imputed_vcf
        File imputed_vcf_index = GlimpseLigate.imputed_vcf_index
        File imputed_vcf_md5sum = GlimpseLigate.imputed_vcf_md5sum
        
        File qc_metrics = CollectQCMetrics.qc_metrics
        File coverage_metrics = CombineCoverageMetrics.coverage_metrics
    }
}

task ComputeShardsAndMemoryPerShard {
    input {
        File reference_chunks_memory
        Array[String] contigs
        Int n_samples
    }

    command <<<
        python3 << EOF
        import pandas as pd
        import numpy as np


        df = pd.read_csv('~{reference_chunks_memory}', sep='\t', header=None, names=['contig', 'reference_shard', 'base_gb', 'slope_per_sample_gb'])

        # filter dataframe by contig list
        chromosomes_to_filter = ["~{sep='", "' contigs}"]
        filtered_df = df[df['contig'].isin(chromosomes_to_filter)]

        # write out reference shards to process
        filtered_df['reference_shard'].to_csv('reference_shard_file_paths.tsv', sep='\t', index=False, header=None)

        # calculate memory usage and save to file
        filtered_df['mem_gb'] = filtered_df['base_gb'] + filtered_df['slope_per_sample_gb'] * ~{n_samples}
        filtered_df['mem_gb'] = filtered_df['mem_gb'].apply(lambda x: min(256, int(np.ceil(x))))  # cap at 256 GB
        filtered_df['mem_gb'].to_csv('memory_per_chunk.tsv', sep='\t', index=False, header=None)
        EOF
    >>>

    runtime {
        docker : "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
    }

    output {
        Array[String] reference_chunk_file_paths = read_lines("reference_shard_file_paths.tsv")
        Array[Int] mem_gb_per_chunk = read_lines("memory_per_chunk.tsv")
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
        File coverage_metrics = "phase_output_stats_coverage.txt.gz"
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
    }

    parameter_meta {
        imputed_vcf: {
                        localization_optional: true
                    }
    }

    Int disk_size_gb = 100
    
    command <<<
        set -euo pipefail

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
