version 1.0

import "../../../tasks/wdl/Utilities.wdl" as utils

workflow MapMyCells {
    meta {
        description: "Executes the MapMyCells cell_type_mapper for cell type annotation."
        allowNestedInputs: true
    }

    String pipeline_version = "0.1.0"

    input {
        File query_h5ad
        String input_id

        # Select the reference atlas. Allowed values: "Human_MTG", "Mouse_WMB", "Custom"
        String reference_atlas = "Human_MTG"

        # If reference_atlas == "Custom", you must provide the following:
        File? custom_precomputed_stats
        File? custom_gene_mapping_db
        File? custom_query_markers

        String algorithm = "hierarchical"
        String normalization = "raw"

        # Docker image
        # ponytail: pinned to the warp-tools add-mapmycells branch tag; repin to a
        # versioned tag once warp-tools/3rd-party-tools/mapmycells/docker_versions.tsv
        # gets a real bump and this pipeline merges.
        String docker = "us.gcr.io/broad-gotc-prod/mapmycells:add-mapmycells"

        # Runtime
        Int cpu = 8
        Int memory_gb = 64
        Int disk_size = 100
        Int preemptible = 3
    }

    if (reference_atlas == "Custom" && !defined(custom_precomputed_stats)) {
        call utils.ErrorWithMessage as ErrorCustomPrecomputedStatsRequired {
            input:
                message = "custom_precomputed_stats is required when reference_atlas is 'Custom'."
        }
    }

    File precomputed_stats = if (reference_atlas == "Human_MTG") then "gs://broad-gotc-test-storage/mapmycells/precomputed_stats_10X_MTG_revision_230821.h5" else if (reference_atlas == "Mouse_WMB") then "gs://broad-gotc-test-storage/mapmycells/precomputed_stats_ABC_revision_230821.h5" else select_first([custom_precomputed_stats])

    File? gene_mapping_db = if (reference_atlas == "Human_MTG" || reference_atlas == "Mouse_WMB") then "gs://broad-gotc-test-storage/mapmycells/mmc_gene_mapper.2025-08-04.db" else custom_gene_mapping_db

    # Only Mouse_WMB ships a canonical markers file today; Human_MTG and Custom both
    # fall back to on-the-fly marker selection unless the caller supplies their own.
    File? query_markers = if (reference_atlas == "Mouse_WMB") then "gs://broad-gotc-test-storage/mapmycells/mouse_markers_230821.json" else custom_query_markers

    call RunMapMyCells {
        input:
            query_h5ad = query_h5ad,
            precomputed_stats = precomputed_stats,
            gene_mapping_db = gene_mapping_db,
            query_markers = query_markers,
            input_id = input_id,
            algorithm = algorithm,
            normalization = normalization,
            docker = docker,
            cpu = cpu,
            memory_gb = memory_gb,
            disk_size = disk_size,
            preemptible = preemptible
    }

    output {
        File output_csv = RunMapMyCells.output_csv
        File output_json = RunMapMyCells.output_json
    }
}

task RunMapMyCells {
    input {
        File query_h5ad
        File precomputed_stats
        File? gene_mapping_db
        File? query_markers
        String input_id
        String algorithm
        String normalization
        String docker
        Int cpu
        Int memory_gb
        Int disk_size
        Int preemptible
    }

    command <<<
        set -euo pipefail
        
        # Determine the entrypoint based on whether markers are provided
        if [ -n "~{query_markers}" ]; then
            python -m cell_type_mapper.cli.from_specified_markers \
                --query_path ~{query_h5ad} \
                --extended_result_path ~{input_id}_mapmycells_extended.json \
                --csv_result_path ~{input_id}_mapmycells_results.csv \
                --type_assignment.n_processors ~{cpu} \
                --precomputed_stats.path ~{precomputed_stats} \
                --query_markers.serialized_lookup ~{query_markers} \
                --query_markers.collapse_markers False \
                --type_assignment.algorithm ~{algorithm} \
                --type_assignment.normalization ~{normalization} \
                ~{"--gene_mapping.db_path " + gene_mapping_db}
        else
            python -m cell_type_mapper.cli.map_to_on_the_fly_markers \
                --query_path ~{query_h5ad} \
                --extended_result_path ~{input_id}_mapmycells_extended.json \
                --csv_result_path ~{input_id}_mapmycells_results.csv \
                --n_processors ~{cpu} \
                --precomputed_stats.path ~{precomputed_stats} \
                --type_assignment.algorithm ~{algorithm} \
                --type_assignment.normalization ~{normalization} \
                ~{"--gene_mapping.db_path " + gene_mapping_db} \
                --query_markers.n_per_utility 15 \
                --reference_markers.log2_fold_min_th 0.5
        fi
    >>>

    output {
        File output_csv = "~{input_id}_mapmycells_results.csv"
        File output_json = "~{input_id}_mapmycells_extended.json"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_size} HDD"
        preemptible: preemptible
    }
    
    parameter_meta {
        query_h5ad: "Input AnnData file containing raw transcriptomics counts."
        precomputed_stats: "HDF5 file containing the precomputed hierarchical taxonomy stats."
        gene_mapping_db: "Optional SQLite database for translating gene symbols to Ensembl IDs."
        query_markers: "Optional JSON file containing predefined marker genes."
        input_id: "Prefix for output files."
        algorithm: "Type assignment algorithm. Options: 'hierarchical', 'hann'. Default: 'hierarchical'."
        normalization: "Normalization method to use for mapping. Options: 'raw', 'log2CPM'. Default: 'raw'."
        preemptible: "Number of times to attempt to run on a preemptible VM."
    }
}
