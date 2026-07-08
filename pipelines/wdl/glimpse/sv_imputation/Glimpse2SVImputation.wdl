version 1.0

import "./MultilevelHierarchicallyPasteVcfsStreaming.wdl" as MultilevelHierarchicallyPasteVcfsStreaming
import "./Glimpse2SVImputationBatched.wdl" as Glimpse2SVImputationBatched

workflow Glimpse2SVImputation {
    String pipeline_version = "0.0.1"
    String batch_pipeline_version = "0.0.1"

    input {
        File? input_gvcfs_fofn
        File? input_gvcf_idxs_fofn
        File? sample_names_file          # order of sample names must match that of gVCFs

        Array[File]? input_gvcfs
        Array[File]? input_gvcf_idxs
        Array[String]? entity_ids
        File? sample_names_map_file           # TSV map of entity_id (research_id) to id2 for AoU DRAGEN gVCFs;
                                              # Terra struggles with id2 as they are parsed as mixed strings/numbers

        String output_prefix

        # inputs for PreprocessPLs
        File preprocess_panel_bubble_split_sites_only_vcf       # can be subset of panel, e.g., simple bubble alleles only
        File preprocess_panel_bubble_split_sites_only_vcf_idx
        File? extract_bubble_likelihoods_script
        File? extract_bubble_likelihoods_cargo_toml
        File? extract_bubble_likelihoods_binary
        String? extract_bubble_likelihoods_extra_args

        File paste_vcfs_binary
        Array[String] paste_regions

        # inputs for Batched wdl
        # inputs for PopAndMarginalizeCollisions
        File pop_glimpse2_panel_resources_json
        File? pop_glimpse2_script               # heavily modified version of convert-to-biallelic.py
        File? pop_glimpse2_cargo_toml
        File? pop_glimpse2_binary

        String chromosome

        File genetic_maps_tsv
        File chunked_panel_json

        File? remap_sample_names_file    # TSV with old_name new_name mappings
    }

    if (defined(input_gvcfs_fofn)) {
        Array[File] parsed_gvcfs = read_lines(select_first([input_gvcfs_fofn]))
    }
    Array[File] input_gvcfs_ = select_first([input_gvcfs, parsed_gvcfs])

    if (defined(input_gvcf_idxs_fofn)) {
        Array[File] parsed_gvcf_idxs = read_lines(select_first([input_gvcf_idxs_fofn]))
    }
    Array[File] input_gvcf_idxs_ = select_first([input_gvcf_idxs, parsed_gvcf_idxs])

    if (defined(sample_names_file)) {
        Array[String] parsed_sample_names = read_lines(select_first([sample_names_file]))
    }

    # Replaced map scatter with a bash task call
    if (defined(entity_ids) && defined(sample_names_map_file)) {
        call MapSampleNames {
            input:
                entity_ids = select_first([entity_ids]),
                sample_names_map_file = select_first([sample_names_map_file])
        }
    }

    Array[String] sample_names_ = select_first([MapSampleNames.mapped_sample_names, parsed_sample_names])

    scatter (j in range(length(input_gvcfs_))) {
        call PreprocessPLs as PreprocessPLsGVCF {
            input:
                input_vcf = input_gvcfs_[j],
                input_vcf_idx = input_gvcf_idxs_[j],
                mode = "gvcf",
                panel_bubble_split_sites_only_vcf = preprocess_panel_bubble_split_sites_only_vcf,
                panel_bubble_split_sites_only_vcf_idx = preprocess_panel_bubble_split_sites_only_vcf_idx,
                sample_names = [sample_names_[j]],
                output_prefix = output_prefix + ".sample-" + j + "." + sample_names_[j] + ".preprocessedPLs",
                extract_bubble_likelihoods_script = extract_bubble_likelihoods_script,
                cargo_toml = extract_bubble_likelihoods_cargo_toml,
                extract_bubble_likelihoods_binary = extract_bubble_likelihoods_binary,
                extra_args = extract_bubble_likelihoods_extra_args
        }
    }

    # two-level localized hierarchical merge over entire chromosome
    call MultilevelHierarchicallyPasteVcfsStreaming.HierarchicallyMergeVcfs as PastePreprocessPLsGVCFs {
        input:
            vcfs_array = PreprocessPLsGVCF.preprocessed_pls_vcf,
            vcf_idxs_array = PreprocessPLsGVCF.preprocessed_pls_vcf_idx,
            regions = paste_regions,
            batch_sizes = [50, 50],
            do_localization = [true, true],
            timeouts_min = [0, 0],
            output_prefix = output_prefix + ".preprocessedPLs",
            paste_vcfs_binary = paste_vcfs_binary,
            extra_merge_args = "--threads $(nproc) --format GT,PL",
            extra_concat_args = "--threads $(nproc) --naive"
    }

    call Glimpse2SVImputationBatched.Glimpse2SVImputationBatched {
        input:
            input_preprocessed_joint_vcf = PastePreprocessPLsGVCFs.merged_vcf,
            input_preprocessed_joint_vcf_idx = PastePreprocessPLsGVCFs.merged_vcf_idx,
            remap_sample_names_file = remap_sample_names_file,
            chromosome = chromosome,
            genetic_maps_tsv = genetic_maps_tsv,
            chunked_panel_json = chunked_panel_json,
            output_prefix = output_prefix,
            pop_glimpse2_panel_resources_json = pop_glimpse2_panel_resources_json,
            pop_glimpse2_script = pop_glimpse2_script,
            pop_glimpse2_cargo_toml = pop_glimpse2_cargo_toml,
            pop_glimpse2_binary = pop_glimpse2_binary
    }

    output {
        File glimpse2_bubble_posteriors_vcf = Glimpse2SVImputationBatched.glimpse2_bubble_posteriors_vcf
        File glimpse2_bubble_posteriors_vcf_idx = Glimpse2SVImputationBatched.glimpse2_bubble_posteriors_vcf_idx
        File glimpse2_popped_posteriors_vcf = Glimpse2SVImputationBatched.glimpse2_popped_posteriors_vcf
        File glimpse2_popped_posteriors_vcf_idx = Glimpse2SVImputationBatched.glimpse2_popped_posteriors_vcf_idx
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

task MapSampleNames {
    input {
        Array[String] entity_ids
        File sample_names_map_file
    }

    command <<<
        set -euo pipefail

        # Use awk to load the TSV map into memory, then translate the entity IDs array in order
        awk 'BEGIN {FS="\t"; OFS="\t"}
        NR==FNR {
            # First pass: read the map file into an array
            map[$1] = $2;
            next
        }
        {
            # Second pass: read the entity_ids file
            if ($1 in map) {
                print map[$1]
            } else {
                print "Error: ID " $1 " not found in map file" > "/dev/stderr"
                exit 1
            }
        }' ~{sample_names_map_file} ~{write_lines(entity_ids)} > mapped_names.txt
    >>>

    output {
        Array[String] mapped_sample_names = read_lines("mapped_names.txt")
    }

    runtime {
        docker: "ubuntu:22.04"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 10 HDD"
        preemptible: 3
    }
}

task PreprocessPLs {
    input {
        File input_vcf
        File input_vcf_idx
        String mode     # joint or gvcf
        File panel_bubble_split_sites_only_vcf
        File panel_bubble_split_sites_only_vcf_idx
        String? output_region
        Array[String] sample_names
        String output_prefix

        File? extract_bubble_likelihoods_script
        File? cargo_toml
        File? extract_bubble_likelihoods_binary
        String? extra_args = "--window 15000 --cap-pl 30 --scale-pl 5.0 --threads $(nproc)"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 2 * ceil(size([input_vcf, panel_bubble_split_sites_only_vcf], "GB"))

    File sample_names_list = write_lines(sample_names)

    command <<<
        set -euxo pipefail

        if [ -n "~{extract_bubble_likelihoods_binary}" ]; then
            EXTRACT_BIN="~{extract_bubble_likelihoods_binary}"
            chmod +x $EXTRACT_BIN
        else
            mkdir -p extract-bubble-PLs/src
            cp ~{extract_bubble_likelihoods_script} extract-bubble-PLs/src/main.rs
            cp ~{cargo_toml} extract-bubble-PLs
            cd extract-bubble-PLs
            cargo build --release
            cd ..
            EXTRACT_BIN="./extract-bubble-PLs/target/release/extract-bubble-PLs"
        fi

        $EXTRACT_BIN ~{mode} \
            ~{panel_bubble_split_sites_only_vcf}##idx##~{panel_bubble_split_sites_only_vcf_idx} \
            ~{input_vcf}##idx##~{input_vcf_idx} \
            ~{output_prefix}.bcf \
            ~{"--region " + output_region} \
            --samples ~{sample_names_list} \
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
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-rust:v1"
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
    }
}
