version 1.0

import "../../../../tasks/wdl/ImputationTasks.wdl" as tasks

workflow QuotaConsumed {
    # if this changes, update the quota_consumed_version value in Glimpse2LowPassImputation.wdl
    String pipeline_version = "0.0.1"

    input {
        Array[String] contigs

        # this is the path to a directory that contains sites vcf, sites table, and reference chunks file. should end with a "/"
        String reference_panel_prefix

        # service currently does not accept VCFs as input
        Array[File]? crams
        Array[File]? cram_indices
        Array[String]? sample_ids
        File? cram_manifest
        File fasta
        File fasta_index
        String output_basename

        File ref_dict
    }

    # validate that either crams, or cram manifest is provided
    if (defined(crams)) {
        Int crams_quota_consumed = length(select_first([crams]))
    }

    if (defined(cram_manifest)) {
        call CountCramsFromManifest {
            input:
                cram_manifest = select_first([cram_manifest])
        }
    }

    output {
        Int quota_consumed = select_first([crams_quota_consumed, CountCramsFromManifest.cram_manifest_count, 0])
    }
}

task CountCramsFromManifest {
    input {
        File cram_manifest

        String docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
        Int cpu = 1
        Int memory_mb = 4000
        Int disk_size_gb = ceil(size(cram_manifest, "GiB")) + 10
    }

    command <<<
        set -e -o pipefail

        grep ".cram" ~{cram_manifest} | wc -l > cram_manifest_count.txt
    >>>

    output {
        Int cram_manifest_count = read_int("cram_manifest_count.txt")
    }
    runtime {
        docker: docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 0
        maxRetries: 1
        noAddress: true
    }
}
