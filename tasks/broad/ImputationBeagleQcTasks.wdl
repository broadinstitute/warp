task QcChecks {
    input {
        File vcf_input

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
        Int cpu = 1
        Int memory_mb = 4000
        Int disk_size_gb = ceil(1.1*size(vcf_input, "GiB")) + 10
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    String vcf_basename = basename(vcf_input)

    command <<<
        # grab header from vcf
        bcftools view -h ~{vcf_input} > header.vcf

        if grep -q "^##fileformat=VCFv[4-9]\." header.vcf; then
            echo "VCF version 4.0+";
        else
            echo "VCF version < 4.0 or not found;" >> qc_messages.txt;
        fi

        gunzip -c ~{vcf_input} | grep -v "#" | cut -f1 | sort -u > chromosomes.txt
        missing_chromosomes=()
        # Check for chr1 through chr22
        for i in {1..22}; do
            if ! grep -q "^chr${i}$" "~{vcf_input}"; then
            missing_chromosomes+=("chr${i}")
            fi
        done
        # report missing chromosomes
        if [ ${#missing_chromosomes[@]} -eq 0 ]; then
            echo "All chromosomes from chr1 to chr22 are present."
        else
            echo "Missing data for chromosomes: ${missing_chromosomes[*]} ;" >> qc_messages.txt
        fi

        # Task should always succeed
        exit 0

    >>>
    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        maxRetries: 2
        noAddress: true
    }
    output {
        String qc_messages = read_string("qc_messages.txt")
    }
}
