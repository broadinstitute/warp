version 1.0

workflow StargazerFromJointVCF {
    input {
        File inputVCF
        File inputVCF_index
        String sample_name
        String outputVCFname = sample_name  + ".pgx_genotyped.vcf.gz"
        Array[String] gene_names

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int? disk_size
    }
    String pipeline_version = "aou_9.0.0"


    call SelectVariants {
        input:
            inputVCF = inputVCF,
            inputVCFindex = inputVCF_index,
            outputVCFname = outputVCFname,
            disk_size = disk_size,
            sample_name = sample_name
    }

    call RunStargazer {
        input:
            input_vcf = SelectVariants.output_vcf,
            input_vcf_index = SelectVariants.output_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            sample_name = sample_name,
            gene_names = gene_names
    }

    call RunStargazer as Stargazer_DPYD {
        input:
            input_vcf = SelectVariants.output_vcf,
            input_vcf_index = SelectVariants.output_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            sample_name = sample_name,
            gene_names = ['dpyd'],
            memory_in_mb = 7500
    }
    
    output {
        Array[File] stargazer_details = flatten([RunStargazer.stargazer_details, Stargazer_DPYD.stargazer_details])
    }
}

    task RunStargazer {
        input {
            File input_vcf
            File input_vcf_index
            File ref_fasta
            File ref_fasta_index
            String sample_name
            Array[String] gene_names
            File? panel_vcf_override
            String? stargazer_docker
            Int? memory_in_mb
        }

    command <<<
        set -euxo pipefail
        ~{if defined(panel_vcf_override) then "mv " + panel_vcf_override + " /stargazer-grc38-v.2.0.2/stargazer/1kgp_vcf/grc38/" else ""}

        # Create REF_CACHE. Used when indexing a CRAM
        /usr/bin/seq_cache_populate.pl -root ./ref/cache ~{ref_fasta} > /dev/null #don't clog the log with 3600 contig names
        export REF_PATH=:
        export REF_CACHE=./ref/cache/%2s/%2s/%s

        gene_names_array=(~{sep=" " gene_names})
        for gene_name in ${gene_names_array[@]}; do

            python /stargazer-grc38-v.2.0.2/stargazer  -t ${gene_name} -o ~{sample_name}_${gene_name} -i ~{input_vcf} -a grc38
    
            mv ~{sample_name}_${gene_name}/report.tsv ~{sample_name}_${gene_name}_report.tsv
            mv ~{sample_name}_${gene_name}/genotype-calls.tsv ~{sample_name}_${gene_name}_genotype-calls.tsv

        done
    >>>

    output {
        Array[File] stargazer_output = glob("~{sample_name}_*_report.tsv")
        Array[File] stargazer_details = glob("~{sample_name}_*_genotype-calls.tsv")
    }

    runtime {
        memory: select_first([memory_in_mb, 3750]) + " MiB"
        disks: "local-disk 100 HDD"
        disk: "100 GB"
        bootDiskSizeGb: 15
        preemptible: 3
        docker: select_first([stargazer_docker, "us.gcr.io/broad-dsde-methods/stargazer:latest"])
    }
}

task SelectVariants {
    input {
        File inputVCF
        File inputVCFindex
        String outputVCFname
        String sample_name

        Int? disk_size = 100
    }

    parameter_meta{
        inputVCF: {localization_optional: true}
        inputVCFindex: {localization_optional: true}
    }

    command <<<
        set -euxo pipefail
        /gatk/gatk SelectVariants -V ~{inputVCF} -sn ~{sample_name} -O ~{outputVCFname}
    >>>
    runtime {
        memory: "7 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        docker: "us.gcr.io/broad-gatk/gatk:4.4.0.0"
    }
    output {
        File output_vcf = "~{outputVCFname}"
        File output_vcf_index = "~{outputVCFname}.tbi"
    }
}