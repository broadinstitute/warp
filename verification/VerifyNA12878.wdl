version 1.0

# Run the Genomes in a Bottle NA12878 VCF evaluation over input vcfs given a validated truth vcf and confidence region
workflow VerifyNA12878 {
    input {
        Array[File] vcf_files   # VCF files to evaluate
        Array[File] vcf_file_indexes
        Array[String] vcf_names

        File truth_vcf = "gs://broad-gotc-test-storage/VerifyNA12878/nist_na12878_giab_hg38_sd_fix.vcf.gz"
        File truth_vcf_index = "gs://broad-gotc-test-storage/VerifyNA12878/nist_na12878_giab_hg38_sd_fix.vcf.gz.tbi"
        File truth_intervals = "gs://broad-gotc-test-storage/VerifyNA12878/HG001_NA12878_GRCh38_GIAB_highconf.exome.interval_list"

        Int? preemptible_attempts
        Int? disk_space
        Int? mem_gb
        Int? cpu
    }

    call RunValidation {
        input:
            vcf_files = vcf_files,
            vcf_file_indexes = vcf_file_indexes,
            vcf_names = vcf_names,
            truth_vcf = truth_vcf,
            truth_vcf_index = truth_vcf_index,
            truth_intervals = truth_intervals,
            cpu = cpu,
            mem_gb = mem_gb,
            disk_space = disk_space,
            preemptible_attempts = preemptible_attempts
    }
    output {
        Array[File] concordance_tables = RunValidation.concordance_tables
    }
}

task RunValidation {
    input {
        Array[File] vcf_files
        Array[File] vcf_file_indexes
        Array[String] vcf_names

        File truth_vcf
        File truth_vcf_index
        File truth_intervals

        # Runtime parameters
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space
        Int? cpu
        Boolean use_ssd = false
    }

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 16000
    Int default_disk_space_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000
    Int max_heap = machine_mem - 500

    #run all the steps here because they're lickety split fast
    command <<<
        files=(~{sep=" " vcf_files})
        names=(~{sep=" " vcf_names})
        for ((i=0;i<${#files[@]};++i)); do
            gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" SelectVariants -V ${files[i]} -sn NA12878 --exclude-non-variants \
            --remove-unused-alternates -O ${names[i]}.NA12878.vcf.gz

            gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" Concordance -eval ${names[i]}.NA12878.vcf.gz \
            --truth ~{truth_vcf} -L ~{truth_intervals} --summary ${names[i]}.summary.tsv
        done
    >>>

    output {
        Array[File] concordance_tables = glob("./*.summary.tsv")
    }

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"

        memory: machine_mem + " MiB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }
}