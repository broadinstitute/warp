version 1.0

task CalculateSomaticContamination {
    input {
        File? intervals
        File reference
        File reference_dict
        File reference_index
        File tumor_cram_or_bam
        File tumor_crai_or_bai
        File? normal_cram_or_bam
        File? normal_crai_or_bai
        File contamination_vcf
        File contamination_vcf_index

        # runtime
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
        File? gatk_override
        Int? additional_disk
        Int memory_mb = 3000
        Int? preemptible_attempts
        Int? max_retries
    }

    Int disk_size = ceil(size(tumor_cram_or_bam,"GB") + size(normal_cram_or_bam,"GB")) + select_first([additional_disk, 10])

    Int command_mem = memory_mb - 500

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx~{command_mem}m" GetPileupSummaries \
             -R ~{reference} \
             -I ~{tumor_cram_or_bam} \
             ~{"--interval-set-rule INTERSECTION -L " + intervals} \
             -V ~{contamination_vcf} \
             -L ~{contamination_vcf} \
             -O pileups.table

        if [[ -f "~{normal_cram_or_bam}" ]];
        then
            gatk --java-options "-Xmx~{command_mem}m" GetPileupSummaries \
                 -R ~{reference} \
                 -I ~{normal_cram_or_bam} \
                 ~{"--interval-set-rule INTERSECTION -L " + intervals} \
                 -V ~{contamination_vcf} \
                 -L ~{contamination_vcf} \
                 -O normal_pileups.table

            gatk --java-options "-Xmx~{command_mem}m" CalculateContamination \
                 -I pileups.table \
                 -O contamination.table \
                 --tumor-segmentation segments.table \
                 -matched normal_pileups.table
        else
            touch normal_pileups.table
            gatk --java-options "-Xmx~{command_mem}m" CalculateContamination \
                 -I pileups.table \
                 -O contamination.table \
                 --tumor-segmentation segments.table
        fi

        grep -v ^sample contamination.table | awk '{print($2)}' > contam.txt
    >>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: memory_mb + " MiB"
        maxRetries: select_first([max_retries, 2])
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
    }

    output {
        File pileups = "pileups.table"
        File normal_pileups = "normal_pileups.table"
        File contamination_table = "contamination.table"
        File maf_segments = "segments.table"
        File contamination = "contam.txt"
    }
}
