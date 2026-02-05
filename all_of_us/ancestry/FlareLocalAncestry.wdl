version 1.0

workflow FlareLocalAncestry {

    String pipeline_version = "1.0.0"

    input {
        File reference_vcf
        File reference_vcf_index
        File plink_map_file
        File ref_panel_mapping_file
        File analysis_vcf
        String output_basename
        File? model

        String flare_docker = "us.gcr.io/broad-gotc-prod/flare:0.0.1"
    }

    call Flare {
        input:
            reference_vcf = reference_vcf,
            reference_vcf_index = reference_vcf_index,
            plink_map_file = plink_map_file,
            ref_panel_mapping_file = ref_panel_mapping_file,
            analysis_vcf = analysis_vcf,
            basename = output_basename,
            model = model,
            flare_docker = flare_docker
    }

    output {
        File log_file = Flare.log_file
        File? model_file = Flare.model_file
        File local_ancestry_vcf = Flare.local_ancestry_vcf
        File global_ancestry_file = Flare.global_ancestry_file
    }

    meta {
        allowNestedInputs: true
    }
}

task Flare {
    input {
        File reference_vcf
        File reference_vcf_index
        File plink_map_file
        File ref_panel_mapping_file
        File analysis_vcf
        String basename

        File? model
        Int seed = 12345

        Int cpu = 1
        Int memory_mb = 6000
        Int disk_size_gb = ceil(size(reference_vcf, "GiB") + size(plink_map_file, "GiB") + size(ref_panel_mapping_file, "GiB") + 2*size(analysis_vcf, "GiB") ) + 10
        String flare_docker
    }

    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    command <<<
        set -euo pipefail

        java -Xms~{command_mem}m -Xmx~{max_heap}m -jar /opt/flare/flare.jar \
        ref=~{reference_vcf} \
        ref-panel=~{ref_panel_mapping_file} \
        map=~{plink_map_file} \
        gt=~{analysis_vcf} \
        seed=~{seed} \
        nthreads=~{cpu} \
        ~{"em=false model=" + model} \
        out=~{basename}

    >>>

    output {
        File log_file = "~{basename}.log"
        File? model_file = "~{basename}.model"
        File local_ancestry_vcf = "~{basename}.anc.vcf.gz)"
        File global_ancestry_file = "~{basename}.global.anc.gz"
    }

    runtime {
        docker: flare_docker
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 0
        maxRetries: 1
        noAddress: true
    }
}
