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

    call RunMapMyCells {
        input:
            query_h5ad               = query_h5ad,
            input_id                 = input_id,
            reference_atlas          = reference_atlas,
            custom_precomputed_stats = custom_precomputed_stats,
            custom_gene_mapping_db   = custom_gene_mapping_db,
            custom_query_markers     = custom_query_markers,
            algorithm                = algorithm,
            normalization            = normalization,
            docker                   = docker,
            cpu                      = cpu,
            memory_gb                = memory_gb,
            disk_size                = disk_size,
            preemptible              = preemptible
    }

    output {
        File output_csv = RunMapMyCells.output_csv
        File output_json = RunMapMyCells.output_json
    }
}

task RunMapMyCells {
    input {
        File query_h5ad
        String input_id
        String reference_atlas
        File? custom_precomputed_stats
        File? custom_gene_mapping_db
        File? custom_query_markers
        String algorithm
        String normalization
        String docker
        Int cpu
        Int memory_gb
        Int disk_size
        Int preemptible
    }

    # Human_MTG and Mouse_WMB reference-atlas assets are baked into the docker image
    # (see warp-tools/3rd-party-tools/mapmycells) instead of hosted externally, so they
    # need no File input / localization of their own -- only their in-container path.
    command <<<
        set -euo pipefail

        case "~{reference_atlas}" in
            Human_MTG)
                precomputed_stats="/opt/mapmycells/data/precomputed_stats.20231120.sea_ad.MTG.h5"
                query_markers=""
                ;;
            Mouse_WMB)
                precomputed_stats="/opt/mapmycells/data/precomputed_stats_ABC_revision_230821.h5"
                query_markers="/opt/mapmycells/data/mouse_markers_230821.json"
                ;;
            Custom)
                precomputed_stats="~{custom_precomputed_stats}"
                query_markers="~{custom_query_markers}"
                ;;
            *)
                >&2 echo "Error: reference_atlas must be one of Human_MTG, Mouse_WMB, Custom (got '~{reference_atlas}')"
                exit 1
                ;;
        esac

        # Determine the entrypoint based on whether markers are available
        if [ -n "$query_markers" ]; then
            python -m cell_type_mapper.cli.from_specified_markers \
                --query_path ~{query_h5ad} \
                --extended_result_path ~{input_id}_mapmycells_extended.json \
                --csv_result_path ~{input_id}_mapmycells_results.csv \
                --type_assignment.n_processors ~{cpu} \
                --precomputed_stats.path "$precomputed_stats" \
                --query_markers.serialized_lookup "$query_markers" \
                --query_markers.collapse_markers False \
                --type_assignment.algorithm ~{algorithm} \
                --type_assignment.normalization ~{normalization} \
                ~{"--gene_mapping.db_path " + custom_gene_mapping_db}
        else
            python -m cell_type_mapper.cli.map_to_on_the_fly_markers \
                --query_path ~{query_h5ad} \
                --extended_result_path ~{input_id}_mapmycells_extended.json \
                --csv_result_path ~{input_id}_mapmycells_results.csv \
                --n_processors ~{cpu} \
                --precomputed_stats.path "$precomputed_stats" \
                --type_assignment.algorithm ~{algorithm} \
                --type_assignment.normalization ~{normalization} \
                ~{"--gene_mapping.db_path " + custom_gene_mapping_db} \
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
        input_id: "Prefix for output files."
        reference_atlas: "Which reference atlas's baked-in assets to use, or 'Custom' to supply your own."
        custom_precomputed_stats: "HDF5 file containing the precomputed hierarchical taxonomy stats. Required when reference_atlas is 'Custom'."
        custom_gene_mapping_db: "Optional SQLite database for translating gene symbols to Ensembl IDs."
        custom_query_markers: "Optional JSON file containing predefined marker genes. Only used when reference_atlas is 'Custom'."
        algorithm: "Type assignment algorithm. Options: 'hierarchical', 'hann'. Default: 'hierarchical'."
        normalization: "Normalization method to use for mapping. Options: 'raw', 'log2CPM'. Default: 'raw'."
        preemptible: "Number of times to attempt to run on a preemptible VM."
    }
}
