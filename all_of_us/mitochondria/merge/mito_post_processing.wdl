version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow MitoPostProcessing {
    meta {
        description: "Runs mito post-processing: exports filtered VCF and sample metadata TSV."
        allowNestedInputs: true
    }

    input {
        String output_path
        String input_path
        String output_base

        String hail_docker = "us.gcr.io/broad-gotc-prod/aou_mitochondria_post:0.0.8"
        RuntimeAttr? runtime_attr_override
    }

    String pipeline_version = "aou_9.0.1"

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
        max_retries:      1
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    command <<<
        set -euo pipefail

        # Hail's Spark/JVM backend spills temporary data to TMPDIR.  In Cromwell,
        # /tmp is on the small OS boot disk, not the provisioned SSD.  Redirect
        # both the JVM tmpdir and Hail's own tmp_dir to the execution directory,
        # which IS on the provisioned disk.
        mkdir -p hail_tmp
        export TMPDIR="$(pwd)/hail_tmp"
        export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$(pwd)/hail_tmp"

        python3 /opt/mito_plot_filter.py \
            --input-path  "~{input_path}" \
            --output-root "~{output_path}" \
            --basename    "~{output_base}" \
            --tmp-dir     "$(pwd)/hail_tmp"
    >>>

    output {
        File filtered_vcf                       = "~{output_base}.filtered.vcf.bgz"
        File filtered_vcf_tbi                   = "~{output_base}.filtered.vcf.bgz.tbi"
        File sample_metadata_tsv                = "~{output_base}.metadata.tsv"
    }

    runtime {
        memory:        select_first([runtime_override.mem_gb,           runtime_default.mem_gb])           + " GB"
        disks:         "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb])   + " SSD"
        cpu:           select_first([runtime_override.cpu_cores,        runtime_default.cpu_cores])
        preemptible:   select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:    select_first([runtime_override.max_retries,      runtime_default.max_retries])
        docker:        hail_docker
    }
}
