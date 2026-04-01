version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow MitoPostProcessing {
    meta {
        description: "Runs mito post-processing from the cleaned notebook: exports filtered VCF, sample metadata TSV, and all generated plots as SVG."
        allowNestedInputs: true
    }

    input {
        String output_path
        String input_path
        String output_base

        String hail_docker = "us.gcr.io/broad-gotc-prod/aou_mitochondria_post:0.0.1"
        RuntimeAttr? runtime_attr_override
    }

    String pipeline_version = "aou_9.0.0"

    call RunMitoPostProcessing {
        input:
            output_path              = output_path,
            input_path               = input_path,
            output_base              = output_base,
            hail_docker              = hail_docker,
            runtime_attr_override    = runtime_attr_override
    }

    output {
        File filtered_vcf                      = RunMitoPostProcessing.filtered_vcf
        File filtered_vcf_tbi                  = RunMitoPostProcessing.filtered_vcf_tbi
        File sample_metadata_tsv               = RunMitoPostProcessing.sample_metadata_tsv

        File variants_per_sample_svg           = RunMitoPostProcessing.variants_per_sample_svg
        File mito_cn_distribution_svg          = RunMitoPostProcessing.mito_cn_distribution_svg
        File variant_allele_frequency_svg      = RunMitoPostProcessing.variant_allele_frequency_svg
        File variant_af_and_allele_fraction_svg = RunMitoPostProcessing.variant_af_and_allele_fraction_svg
        File numt_fp_by_mtcn_svg               = RunMitoPostProcessing.numt_fp_by_mtcn_svg
        File haplogroup_heteroplasmy_svg       = RunMitoPostProcessing.haplogroup_heteroplasmy_svg
        File haplogroup_homoplasmy_svg         = RunMitoPostProcessing.haplogroup_homoplasmy_svg
    }
}

task RunMitoPostProcessing {
    input {
        String output_path
        String input_path
        String output_base

        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb:           32,
        disk_gb:          200,
        cpu_cores:        8,
        preemptible_tries: 0,
        max_retries:      1,
        boot_disk_gb:     20
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    command <<<
        set -euo pipefail

        python3 /opt/mito_plot_filter.py \
            --input-path  "~{input_path}" \
            --output-path "~{output_path}" \
            --output-base "~{output_base}"
    >>>

    output {
        File filtered_vcf                       = "~{output_base}.vcf.bgz"
        File filtered_vcf_tbi                   = "~{output_base}.vcf.bgz.tbi"
        File sample_metadata_tsv                = "~{output_base}_metadata.tsv"

        File variants_per_sample_svg            = "~{output_base}.variants_per_sample.svg"
        File mito_cn_distribution_svg           = "~{output_base}.mito_cn_distribution.svg"
        File variant_allele_frequency_svg       = "~{output_base}.variant_allele_frequency.svg"
        File variant_af_and_allele_fraction_svg = "~{output_base}.variant_af_and_allele_fraction.svg"
        File numt_fp_by_mtcn_svg                = "~{output_base}.numt_fp_by_mtcn.svg"
        File haplogroup_heteroplasmy_svg        = "~{output_base}.haplogroup_heteroplasmy.svg"
        File haplogroup_homoplasmy_svg          = "~{output_base}.haplogroup_homoplasmy.svg"
    }

    runtime {
        memory:        select_first([runtime_override.mem_gb,           runtime_default.mem_gb])           + " GB"
        disks:         "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb])   + " SSD"
        cpu:           select_first([runtime_override.cpu_cores,        runtime_default.cpu_cores])
        preemptible:   select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:    select_first([runtime_override.max_retries,      runtime_default.max_retries])
        docker:        hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb,   runtime_default.boot_disk_gb])
    }
}
