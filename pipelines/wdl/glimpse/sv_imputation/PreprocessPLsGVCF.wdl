version 1.0

import "./MultilevelHierarchicallyPasteVcfsStreaming.wdl" as MultilevelHierarchicallyPasteVcfsStreaming
import "../../../../tasks/wdl/Glimpse2SVImputationTasks.wdl" as Glimpse2SVImputationTasks

workflow PreprocessPLsGVCF {
    # if this changes, update the preprocessing_pls_gvcf_pipeline_version value in Glimpse2SVImputation.wdl
    String pipeline_version = "0.0.13"
    String multi_level_paste_pipeline_version = "0.0.8"
    input {
        File input_gvcf_manifest

        # inputs for PreprocessPLs
        File preprocess_panel_bubble_split_sites_only_vcf       # can be subset of panel, e.g., simple bubble alleles only
        File preprocess_panel_bubble_split_sites_only_vcf_idx
        String? extract_bubble_likelihoods_extra_args

        Array[String] paste_regions
    }

    call Glimpse2SVImputationTasks.ParseVcfManifestIntoArrays as ParseInputManifest {
        input:
            gvcf_manifest = input_gvcf_manifest
    }

    scatter (j in range(length(ParseInputManifest.input_gvcfs))) {
        call PreprocessPLs as PreprocessPLsGVCF {
            input:
                input_vcf = ParseInputManifest.input_gvcfs[j],
                input_vcf_idx = ParseInputManifest.input_gvcf_idxs[j],
                mode = "gvcf",
                panel_bubble_split_sites_only_vcf = preprocess_panel_bubble_split_sites_only_vcf,
                panel_bubble_split_sites_only_vcf_idx = preprocess_panel_bubble_split_sites_only_vcf_idx,
                output_prefix = "sample-" + j,
                extra_args = extract_bubble_likelihoods_extra_args
        }
    }

    # two-level localized hierarchical merge over entire chromosome
    call MultilevelHierarchicallyPasteVcfsStreaming.MultilevelHierarchicallyMergeVcfs as PastePreprocessPLsGVCFs {
        input:
            vcfs_array = PreprocessPLsGVCF.preprocessed_pls_vcf,
            vcf_idxs_array = PreprocessPLsGVCF.preprocessed_pls_vcf_idx,
            regions = paste_regions,
            batch_sizes = [50, 50],
            do_localization = [true, true],
            timeouts_min = [0, 0],
            output_prefix = "preprocessedPLs.merged",
            extra_merge_args = "--format GT,PL",
            extra_concat_args = "--naive"
    }

    output {
        File preprocessed_pls_vcf = PastePreprocessPLsGVCFs.merged_vcf
        File preprocessed_pls_vcf_idx = PastePreprocessPLsGVCFs.merged_vcf_idx
        Int num_samples = length(ParseInputManifest.input_gvcfs)
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Boolean? use_ssd
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

task PreprocessPLs {
    input {
        File input_vcf
        File input_vcf_idx
        String mode     # joint or gvcf
        File panel_bubble_split_sites_only_vcf
        File panel_bubble_split_sites_only_vcf_idx
        String? output_region
        String output_prefix

        String? extra_args = "--window 15000 --cap-pl 30 --scale-pl 5.0"
        Int cpu = 1

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = ceil(2*size([input_vcf, panel_bubble_split_sites_only_vcf], "GB")) + 10

    command <<<
        set -euxo pipefail

        # Extract sample name from GVCF header and write to file for extract-bubble-PLs tool
        bcftools query -l ~{input_vcf} > sample_name.txt

        /usr/local/bin/extract-bubble-PLs ~{mode} \
            ~{panel_bubble_split_sites_only_vcf}##idx##~{panel_bubble_split_sites_only_vcf_idx} \
            ~{input_vcf}##idx##~{input_vcf_idx} \
            ~{output_prefix}.bcf \
            ~{"--region " + output_region} \
            --samples sample_name.txt \
            --threads ~{cpu} \
            ~{extra_args}

        bcftools index ~{output_prefix}.bcf

        echo "Number of bubble alleles extracted..."
        bcftools index -n ~{output_prefix}.bcf
    >>>

    output {
        File preprocessed_pls_vcf = "~{output_prefix}.bcf"
        File preprocessed_pls_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpu,
        mem_gb:             2,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       0,
        use_ssd:            true,
        preemptible_tries:  4,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/sv-imputation-rust-tools:1.0.0-5dc0f19-1784328222"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        noAddress: true
    }
}
