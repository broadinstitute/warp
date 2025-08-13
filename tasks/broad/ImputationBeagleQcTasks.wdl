version 1.0

task QcChecks {
    input {
        File vcf_input
        File ref_dict

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
        Int cpu = 1
        Int memory_mb = 4000
        Int disk_size_gb = ceil(1.1*size(vcf_input, "GiB")) + 10
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    String vcf_basename = basename(vcf_input)

    command <<<
        # create empty qc messages file
        touch qc_messages.txt

        # index vcf - necessary for gatk command
        bcftools index -t ~{vcf_input}

        # grab header from vcf
        bcftools view -h ~{vcf_input} > header.vcf

        # check for vcf version > 4.0
        if grep -q "^##fileformat=VCFv[4-9]\." header.vcf; then
            echo "VCF version 4.0+";
        else
            echo "VCF version < 4.0 or not found;" >> qc_messages.txt;
        fi

        # check for variants in all canonical chromosomes - chr1 to chr22
        gunzip -c ~{vcf_input} | grep -v "#" | cut -f1 | sort -u > chromosomes.txt

        missing_chromosomes=()
        # Check for chr1 through chr22
        for i in {1..22}; do
            if ! grep -q "^chr${i}$" "chromosomes.txt"; then
                missing_chromosomes+=("chr${i}")
            fi
        done

        if [ ${#missing_chromosomes[@]} -eq 0 ]; then
            echo "All chromosomes from chr1 to chr22 are present."
        else
            echo "Missing data for chromosomes: ${missing_chromosomes[*]};" >> qc_messages.txt
        fi

        # check for a large number of variants in input vcf
        line_count=$(gunzip -c ~{vcf_input} | grep -v "#" | wc -l | tr -d ' ')
        if [ "$line_count" -gt 10000000 ]; then
            echo "Greater than 10 million variants found in input VCF." >> qc_messages.txt
        else
            echo "Less than or equal to 10 million variants found in input VCF."
        fi

        # check reference header lines if they exist
        gatk ValidateVariants \
        -V ~{vcf_input} \
        --sequence-dictionary ~{ref_dict} \
        --validation-type-to-exclude ALL \
        2> gatk_output.txt

        if grep -q "incompatible contigs" gatk_output.txt; then
            echo "found incompatible contigs in VCF header;" >> qc_messages.txt;
        else
            echo "no incompatible contigs found in VCF header;"
        fi

        # This task should always succeed
        exit 0

    >>>
    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 3
        maxRetries: 2
        noAddress: true
    }
    output {
        String qc_messages = read_string("qc_messages.txt")
    }
}
