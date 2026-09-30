version 1.0

task tensorqtl_cis_nominal {

    input {
        File plink_pgen
        File plink_pvar
        File plink_psam

        File phenotype_bed
        File covariates
        String prefix

        File? interaction
        File? phenotype_groups

        Int memory
        Int disk_space
        Int num_threads
        Int num_gpus
        Int num_preempt

        String pipeline_version = "aou_9.0.0"
    }

    command <<<
        set -euo pipefail

        plink_base=$(echo "~{plink_pgen}" | rev | cut -f 2- -d '.' | rev)

        python3 -m tensorqtl \
            $plink_base ~{phenotype_bed} ~{prefix} \
            --mode cis_nominal \
            --covariates ~{covariates} \
            ~{if defined(interaction) then "--interaction " + interaction else ""} \
            ~{if defined(phenotype_groups) then "--phenotype_groups " + phenotype_groups else ""}
    >>>

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/tensorqtl@sha256:f6efb9e592eb32c46cb75070be2769b34381d60cbb2709d2885771324abfe32a"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        bootDiskSizeGb: 25
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
        gpuType: "nvidia-tesla-t4"
        gpuCount: "~{num_gpus}"
        zones: ["us-central1-c"]
    }

    output {
        Array[File] chr_parquet=glob("${prefix}*.parquet")
        File log=glob("${prefix}*.log")[0]
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow tensorqtl_cis_nominal_workflow {
    input {
        File plink_pgen
        File plink_pvar
        File plink_psam
        File phenotype_bed
        File covariates
        String prefix

        File? interaction
        File? phenotype_groups

        Int memory
        Int disk_space
        Int num_threads
        Int num_gpus
        Int num_preempt
    }

    String pipeline_version = "aou_9.0.0"

    call tensorqtl_cis_nominal {
        input:
            plink_pgen = plink_pgen,
            plink_pvar = plink_pvar,
            plink_psam = plink_psam,
            phenotype_bed = phenotype_bed,
            covariates = covariates,
            prefix = prefix,
            interaction = interaction,
            phenotype_groups = phenotype_groups,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            num_gpus = num_gpus,
            num_preempt = num_preempt,
            pipeline_version = pipeline_version
    }

    output {
        Array[File] chr_parquet = tensorqtl_cis_nominal.chr_parquet
        File log = tensorqtl_cis_nominal.log
    }
}
