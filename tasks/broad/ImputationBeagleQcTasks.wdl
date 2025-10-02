version 1.0

task QcChecks {
    input {
        File vcf_input
        Array[String] allowed_contigs
        File ref_dict

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
        Int cpu = 1
        Int memory_mb = 4000
        Int disk_size_gb = ceil(1.1*size(vcf_input, "GiB")) + 10
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    String vcf_basename = basename(vcf_input)

    String ref_dict_basename = basename(ref_dict)

    command <<<
        # create empty qc messages file
        touch qc_messages.txt

        # check for a large number of variants in input vcf and exit if greater than 10 million
        line_count=$(bcftools stats ~{vcf_input}  | grep "number of records:" | awk '{print $6}')
        if [ "$line_count" -gt 10000000 ]; then
            echo "Greater than 10 million variants found in input VCF." >> qc_messages.txt
            exit 0
        else
            echo "Less than or equal to 10 million variants found in input VCF."
        fi

        # grab header from vcf
        bcftools view -h ~{vcf_input} > header.vcf

        # check for vcf version > 4.0
        if grep -q "^##fileformat=VCFv[4-9]\." header.vcf; then
            echo "VCF version 4.0+";
        else
            echo "VCF version < 4.0 or not found." >> qc_messages.txt;
        fi

        # check for variants in at least one of the canonical chromosomes - chr1 to chr22
        gunzip -c ~{vcf_input} | grep -v "#" | cut -f1 | sort -u > chromosomes.txt

        allowed_chromosomes=()
        filtered_chromosomes=()
        # Check for at least one input chromosome in the list of allowed contigs
        for chr in  ~{sep=" " allowed_contigs}; do
            allowed_chromosomes+=("${chr}")
            if grep -q "^${chr}$" "chromosomes.txt"; then
                filtered_chromosomes+=("${chr}")
            fi
        done

        if [ ${#filtered_chromosomes[@]} -eq 0 ]; then
            echo "Input must include data for at least one chromosome in the allowed contigs (${allowed_chromosomes[*]})." >> qc_messages.txt
        else
            echo "Found data for chromosomes: ${filtered_chromosomes[*]}."
        fi

        # check for sorted or non bgzf compressed vcf
        bcftools index -t ~{vcf_input} 2> index_stderr.txt
        if grep -qiE "unsorted positions|not continuous" index_stderr.txt; then
            echo "Input VCF is not sorted" >> qc_messages.txt;
        fi

        if grep -q "not BGZF compressed" index_stderr.txt; then
            echo "Input VCF is BGZF compressed" >> qc_messages.txt;
        fi

        # check reference header lines if they exist
        gatk ValidateVariants \
        -V ~{vcf_input} \
        --sequence-dictionary ~{ref_dict} \
        --validation-type-to-exclude ALL \
        2> gatk_output.txt

        ref_dict_basename="~{ref_dict_basename}"
        if grep -q "incompatible contigs" gatk_output.txt; then
            echo "Found incompatible contigs (against reference dictionary $ref_dict_basename) in VCF header." >> qc_messages.txt;
        else
            echo "No incompatible contigs found in VCF header."
        fi

        # passes_qc is true if qc_messages is empty
        if [ ! -s qc_messages.txt ]; then
            echo "true" > passes_qc.txt
        else
            echo "false" > passes_qc.txt
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
        maxRetries: 1
        noAddress: true
    }
    output {
        Boolean passes_qc = read_boolean("passes_qc.txt")
        String qc_messages = read_string("qc_messages.txt")
    }
}
