version 1.0

workflow QuotaConsumed {
    # if this changes, update the quota_consumed_version value in Glimpse2SVImputation.wdl
    String pipeline_version = "0.0.1"

    input {
        # service expects only gvcf_manifest even though main wdl can alternatively take input arrays
        File gvcf_manifest
        String output_prefix

        # remaining inputs kept for interface consistency with Glimpse2SVImputation.wdl; not all are used by this wdl
        File preprocess_panel_bubble_split_sites_only_vcf
        File preprocess_panel_bubble_split_sites_only_vcf_idx

        Array[String] paste_regions

        Array[String] chromosomes
        File genetic_maps_tsv
        File ref_dict
        File chunked_panel_json

        File pop_glimpse2_panel_resources_json

        Float? info_filter_for_inclusion

        # optional additional header line to add to the output VCF
        String? pipeline_header_line
    }

    call CountGvcfsFromManifest {
        input:
            gvcf_manifest = gvcf_manifest
    }

    output {
        Int quota_consumed = CountGvcfsFromManifest.gvcf_manifest_count
    }
}

task CountGvcfsFromManifest {
    input {
        File gvcf_manifest

        String docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
        Int cpu = 1
        Int memory_mb = 4000
        Int disk_size_gb = ceil(size(gvcf_manifest, "GiB")) + 10
    }

    command <<<
        set -e -o pipefail

        grep -E "\.(vcf|gvcf)\.gz" ~{gvcf_manifest} | wc -l > gvcf_manifest_count.txt
    >>>

    output {
        Int gvcf_manifest_count = read_int("gvcf_manifest_count.txt")
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
