version 1.0

workflow mt_coverage_merge {

    input {
        # Side inputs with fields that are not incoporated in the samples table
        File coverage_tsv
        File ancestry_tsv
        File dob_tsv
        File wgs_median_coverage_tsv 

        File full_data_tsv
        File? sample_list_tsv

        # v2 coverage DB builder parameters
        Boolean skip_summary = false

        # Step 3 (combine vcfs + homref from covdb)
        String vcf_col_name = "final_vcf"
        String combined_mt_name = "combined_vcf"

        # Step 3 sharding controls (v3-style scalable ingestion)
        Boolean shard_step3 = true
        Int step3_shard_size = 2500
        Int step3_merge_fanin = 10
        Int step3_shard_n_partitions = 128

    }

    # (kept for provenance in earlier iterations; currently unused)

    if (defined(sample_list_tsv)) {
        call subset_data_table {
            input:
                full_data_tsv = full_data_tsv,
                sample_list_tsv = sample_list_tsv
            }
    }

    File input_table = select_first([subset_data_table.subset_tsv, full_data_tsv])


    call process_tsv_files {
        input:
            coverage_tsv = coverage_tsv,
            ancestry_tsv = ancestry_tsv,
            dob_tsv = dob_tsv,
            wgs_median_coverage_tsv = wgs_median_coverage_tsv,
            input_tsv = input_table
    }

    call annotate_coverage {
        input:
            input_tsv = process_tsv_files.processed_tsv,  # Input TSV file path
            skip_summary = skip_summary
    }

    if (shard_step3) {
        call make_vcf_shards_from_tsv {
            input:
                input_tsv = process_tsv_files.processed_tsv,
                vcf_col_name = vcf_col_name,
                shard_size = step3_shard_size
        }

        scatter (shard_tsv in make_vcf_shards_from_tsv.shard_tsvs) {
            call build_vcf_shard_mt {
                input:
                    shard_tsv = shard_tsv,
                    n_final_partitions = step3_shard_n_partitions
            }
        }

        # Round 1: merge shard MTs to fewer intermediate MTs (fan-in)
        call make_mt_merge_groups as make_merge_groups_1 {
            input:
                mt_tars = build_vcf_shard_mt.shard_mt_tar,
                fanin = step3_merge_fanin
        }

        scatter (mt_list_tsv in make_merge_groups_1.mt_list_tsvs) {
            call merge_mt_shards as merge_round_1 {
                input:
                    mt_list_tsv = mt_list_tsv,
                    out_mt_name = "merge_round_1_" + basename(mt_list_tsv, ".tsv") + ".mt",
                    chunk_size = step3_merge_fanin
            }
        }

        # Round 2: merge round-1 outputs
        call make_mt_merge_groups as make_merge_groups_2 {
            input:
                mt_tars = merge_round_1.merged_mt_tar,
                fanin = step3_merge_fanin
        }

        Boolean do_merge_round_2 = length(make_merge_groups_2.mt_list_tsvs) > 1

        if (do_merge_round_2) {
            scatter (mt_list_tsv in make_merge_groups_2.mt_list_tsvs) {
                call merge_mt_shards as merge_round_2 {
                    input:
                        mt_list_tsv = mt_list_tsv,
                        out_mt_name = "merge_round_2_" + basename(mt_list_tsv, ".tsv") + ".mt",
                        chunk_size = step3_merge_fanin
                }
            }
        }

        # Round 3: merge round-2 outputs (typically small count, but keep general)
        if (do_merge_round_2) {
            # Round 3 only makes sense if round 2 actually ran.
            call make_mt_merge_groups as make_merge_groups_3 {
                input:
                    mt_tars = select_first([merge_round_2.merged_mt_tar]),
                    fanin = step3_merge_fanin
            }

            Boolean do_merge_round_3 = length(make_merge_groups_3.mt_list_tsvs) > 1

            if (do_merge_round_3) {
                scatter (mt_list_tsv in make_merge_groups_3.mt_list_tsvs) {
                    call merge_mt_shards as merge_round_3 {
                        input:
                            mt_list_tsv = mt_list_tsv,
                            out_mt_name = "merge_round_3_" + basename(mt_list_tsv, ".tsv") + ".mt",
                            chunk_size = step3_merge_fanin
                    }
                }
            }
        }

        # Finalize: apply covdb homref/DP + artifact filter once on the final merged MT
        # Finalize on the deepest merge output that exists.
        if (do_merge_round_2) {
            if (defined(merge_round_3.merged_mt_tar)) {
                call finalize_mt_with_covdb as finalize_mt_with_covdb_round3 {
                    input:
                        in_mt_tar = select_first([merge_round_3.merged_mt_tar])[0],
                        coverage_db_tar = annotate_coverage.output_ht,
                        file_name = combined_mt_name
                }
            }
            if (!defined(merge_round_3.merged_mt_tar)) {
                call finalize_mt_with_covdb as finalize_mt_with_covdb_round2 {
                    input:
                        in_mt_tar = select_first([merge_round_2.merged_mt_tar])[0],
                        coverage_db_tar = annotate_coverage.output_ht,
                        file_name = combined_mt_name
                }
            }
        }

        if (!do_merge_round_2) {
            call finalize_mt_with_covdb as finalize_mt_with_covdb_round1 {
                input:
                    in_mt_tar = merge_round_1.merged_mt_tar[0],
                    coverage_db_tar = annotate_coverage.output_ht,
                    file_name = combined_mt_name
            }
        }
    }

    if (!shard_step3) {
        call combine_vcfs_and_homref_from_covdb {
            input:
                input_tsv = process_tsv_files.processed_tsv,
                coverage_db_tar = annotate_coverage.output_ht,
                vcf_col_name = vcf_col_name,
                file_name = combined_mt_name
        }
    }

    #call add_annotations as annotated {
    #    input:
    #        coverage_mt_tar = annotate_coverage.output_ht,  # Tar.gzipped directory of the Hail table
    #        coverage_tsv = process_tsv_files.processed_tsv,  # Path to the coverage input TSV file
    #        vcf_mt = combine_vcfs.results_tar,  # Path to the MatrixTable
    #        keep_all_samples = true,
    #        output_name = "annotated"
    #}


    #call add_annotations as filt_annotated {
    #    input:
    #        coverage_mt_tar = annotate_coverage.output_ht,  # Tar.gzipped directory of the Hail table
    #        coverage_tsv = process_tsv_files.processed_tsv,  # Path to the coverage input TSV file
    #        vcf_mt = combine_vcfs.results_tar,  # Path to the MatrixTable
    #        keep_all_samples = false,
    #        output_name = "filt_annotated"
    #}

    output {
        File processed_tsv = process_tsv_files.processed_tsv
        File output_coverage_ht = annotate_coverage.output_ht
        File combined_mt_tar = select_first([
            finalize_mt_with_covdb_round3.results_tar,
            finalize_mt_with_covdb_round2.results_tar,
            finalize_mt_with_covdb_round1.results_tar,
            combine_vcfs_and_homref_from_covdb.results_tar
        ])
        #File combined_vcf = combine_vcfs.results_tar
        #File annotated_output_tar = annotated.annotated_output_tar
        #File filt_annotated_output_tar = filt_annotated.annotated_output_tar
    }
}

task make_vcf_shards_from_tsv {
    input {
        File input_tsv
        String vcf_col_name = "final_vcf"
        Int shard_size = 2500

        # Runtime parameters
        Int memory_gb = 16
        Int cpu = 2
        Int disk_gb = 50
    }

    command <<<
        set -euxo pipefail

        # Write shard TSVs directly into the task root.
        # This is the most reliable pattern for Cromwell delocalization across backends.

        python3 <<'EOF'
        import math
        from pathlib import Path

        import pandas as pd

        input_tsv = "~{input_tsv}"
        sample_id_col = "s"
        vcf_col = "~{vcf_col_name}"
        shard_size = int("~{shard_size}")
        out_dir = Path(".")

        if shard_size <= 0:
            raise ValueError("shard_size must be > 0")

        df = pd.read_csv(input_tsv, sep="\t", dtype=str)

        if sample_id_col not in df.columns:
            raise ValueError(f"Missing sample column '{sample_id_col}'")
        if vcf_col not in df.columns:
            raise ValueError(f"Missing VCF column '{vcf_col}'")

        pairs = (
            df[[sample_id_col, vcf_col]]
            .rename(columns={sample_id_col: "s", vcf_col: "vcf"})
            .dropna()
        )
        pairs = pairs[pairs["vcf"].astype(str).str.len() > 0]
        pairs = pairs.sort_values("s", kind="mergesort")

        n_samples = int(pairs.shape[0])
        if n_samples == 0:
            raise ValueError(
                "No usable (s, vcf) rows found. Check that the TSV has non-empty VCF paths in the requested vcf_col."
            )

        n_shards = int(math.ceil(n_samples / shard_size))

        out_dir.mkdir(parents=True, exist_ok=True)
        shard_paths = []
        for shard_idx in range(n_shards):
            start = shard_idx * shard_size
            end = min((shard_idx + 1) * shard_size, n_samples)
            shard_df = pairs.iloc[start:end]
            shard_path = out_dir / f"vcf_shard_{shard_idx:05d}.tsv"
            shard_df.to_csv(shard_path, sep="\t", index=False)
            shard_paths.append(str(shard_path))

        idx_path = out_dir / "shards.tsv"
        with idx_path.open("w") as out:
            for p in shard_paths:
                out.write(f"{p}\n")

        print(f"Wrote {len(shard_paths)} shards")
        EOF

        # Declare shard TSV outputs.
        ls -1 vcf_shard_*.tsv | sort > shard_tsvs.list
    >>>

    output {
        Array[File] shard_tsvs = glob("vcf_shard_*.tsv")
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/python-numpy-pandas:1.0.0-2.2.3-1.25.2"
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " HDD"
    }
}

task build_vcf_shard_mt {
    input {
        File shard_tsv
        Int chunk_size = 100
        Int n_final_partitions = 128
        Boolean overwrite = false
        Boolean include_extra_v2_fields = true

        # Runtime parameters
        Int memory_gb = 64
        Int cpu = 16
        Int disk_gb = 500
        String disk_type = "SSD"
    }

    String out_mt_dirname = basename(shard_tsv, ".tsv") + ".mt"

    # Unique tarball name per shard to avoid collisions.
    String out_tar_name = basename(shard_tsv, ".tsv") + "_shard_mt.tar.gz"

    command <<<
        set -euxo pipefail

        mkdir -p ./tmp
        mkdir -p ./results

        # Spark/Hail scratch
        export SPARK_LOCAL_DIRS="$PWD/tmp"

        DRIVER_MEM_GB=$((~{memory_gb} - 8))
        if [ "$DRIVER_MEM_GB" -lt 4 ]; then DRIVER_MEM_GB=4; fi

        export SPARK_DRIVER_MEMORY="${DRIVER_MEM_GB}g"
        export PYSPARK_SUBMIT_ARGS="--driver-memory ${DRIVER_MEM_GB}g --executor-memory ${DRIVER_MEM_GB}g pyspark-shell"
        export JAVA_OPTS="-Xms${DRIVER_MEM_GB}g -Xmx${DRIVER_MEM_GB}g"

        python3 /opt/mtSwirl/generate_mtdna_call_mt/Terra/build_vcf_shard_mt.py \
            --shard-tsv ~{shard_tsv} \
            --out-mt ./results/~{out_mt_dirname} \
            --temp-dir ./tmp \
            --chunk-size ~{chunk_size} \
            --n-final-partitions ~{n_final_partitions} \
            ~{if overwrite then "--overwrite" else ""} \
            ~{if include_extra_v2_fields then "--include-extra-v2-fields" else ""}

        # Pack as a tar for WDL artifact portability (unique per shard)
        tar -czf "~{out_tar_name}" -C ./results "~{out_mt_dirname}"
    >>>

    output {
        File shard_mt_tar = out_tar_name
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/aou-mitochondrial-combine-vcfs-covdb:dev-hailkeyfix"
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " " + disk_type
    }
}

task make_mt_merge_groups {
    input {
        Array[File] mt_tars
        Int fanin = 10

        # Runtime parameters
        Int memory_gb = 8
        Int cpu = 1
        Int disk_gb = 20
    }

    command <<<
        set -euxo pipefail

        # Serialize the Array[File] into a newline-delimited string for Python.
        # IMPORTANT: export it so the heredoc Python process inherits it.
        # Using newlines avoids shell word-splitting surprises.
        export MT_TARS=$'~{sep="\\n" mt_tars}'
        echo "$MT_TARS"

        python3 <<'EOF'
        import math
        import os

        mt_tars = os.environ.get("MT_TARS", "").strip().splitlines()
        mt_tars = [p for p in mt_tars if p]
        fanin = int("~{fanin}")

        if fanin <= 0:
            raise ValueError("fanin must be > 0")

        n = len(mt_tars)
        if n == 0:
            raise ValueError("mt_tars is empty")

        n_groups = int(math.ceil(n / fanin))
        for g in range(n_groups):
            chunk = mt_tars[g * fanin : (g + 1) * fanin]
            list_path = f"mt_list_{g:05d}.tsv"
            with open(list_path, "w") as out:
                out.write("mt_tar\n")
                for p in chunk:
                    out.write(p + "\n")
        EOF
    >>>

    output {
        Array[File] mt_list_tsvs = glob("mt_list_*.tsv")
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/python-numpy-pandas:1.0.0-2.2.3-1.25.2"
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " HDD"
    }
}

task merge_mt_shards {
    input {
        File mt_list_tsv
        String out_mt_name
        Int chunk_size = 10
        Int min_partitions = 256
        Boolean overwrite = false

        # Runtime parameters
        Int memory_gb = 128
        Int cpu = 16
        Int disk_gb = 1500
        String disk_type = "SSD"
    }

    command <<<
        set -euxo pipefail

        mkdir -p ./tmp
        mkdir -p ./results

        export SPARK_LOCAL_DIRS="$PWD/tmp"

        DRIVER_MEM_GB=$((~{memory_gb} - 8))
        if [ "$DRIVER_MEM_GB" -lt 4 ]; then DRIVER_MEM_GB=4; fi

        export SPARK_DRIVER_MEMORY="${DRIVER_MEM_GB}g"
        export PYSPARK_SUBMIT_ARGS="--driver-memory ${DRIVER_MEM_GB}g --executor-memory ${DRIVER_MEM_GB}g pyspark-shell"
        export JAVA_OPTS="-Xms${DRIVER_MEM_GB}g -Xmx${DRIVER_MEM_GB}g"

        # Localize and extract all input MT tars.
        # NOTE: mt_list_tsv may contain gs:// URIs. We download them *inside the task* (not as inputs)
        # to avoid localizing all tarballs for every task shard.
        mkdir -p ./inputs

        # Sanity-check that required cloud tooling exists (baked into our :dev image)
        command -v gcloud

        # mt_list_tsv has a header: mt_tar
        # For each tar path (usually gs://...), download locally, untar, and write a manifest of local *.mt dirs.
        printf "mt_path\n" > ./inputs/mt_paths.tsv

        i=0
        tail -n +2 ~{mt_list_tsv} | while IFS=$'\t' read -r mt_tar extra; do
          if [ -z "${mt_tar}" ]; then
            continue
          fi

          local_tar="./inputs/mt_$(printf '%05d' ${i}).tar.gz"
          dest_dir="./inputs/mt_$(printf '%05d' ${i}).extract"
          mkdir -p "${dest_dir}"

          if [[ "${mt_tar}" == gs://* ]]; then
            gcloud storage cp "${mt_tar}" "${local_tar}"
          else
            cp -f "${mt_tar}" "${local_tar}"
          fi

          tar -xzf "${local_tar}" -C "${dest_dir}"

          mt_dir=$(find "${dest_dir}" -maxdepth 2 -type d -name "*.mt" | head -n 1)
          if [ -z "${mt_dir}" ]; then
            echo "ERROR: could not find .mt directory after extracting ${local_tar}" >&2
            find "${dest_dir}" -maxdepth 3 -type d | head -100 >&2
            exit 1
          fi
          printf "%s\n" "${mt_dir}" >> ./inputs/mt_paths.tsv

          i=$((i+1))
        done

        python3 /opt/mtSwirl/generate_mtdna_call_mt/Terra/merge_mt_shards.py \
            --mt-list-tsv ./inputs/mt_paths.tsv \
            --out-mt ./results/~{out_mt_name} \
            --temp-dir ./tmp \
            --chunk-size ~{chunk_size} \
            --min-partitions ~{min_partitions} \
            ~{if overwrite then "--overwrite" else ""}

        # Use a unique tarball name per call to avoid collisions across scatter shards and merge rounds.
        tar -czf "~{out_mt_name}.tar.gz" -C ./results "~{out_mt_name}"
    >>>

    output {
        File merged_mt_tar = "~{out_mt_name}.tar.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/aou-mitochondrial-combine-vcfs-covdb:dev-hailkeyfix"
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " " + disk_type
    }
}

task finalize_mt_with_covdb {
    input {
        File in_mt_tar
        File coverage_db_tar
        String artifact_prone_sites_path = "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed"
        String artifact_prone_sites_reference = "default"
        String file_name
        Int minimum_homref_coverage = 100
        Int homref_position_block_size = 1024
        Int n_final_partitions = 1000
        Boolean overwrite = false

        # Runtime parameters
        Int memory_gb = 256
        Int cpu = 32
        Int disk_gb = 2000
        String disk_type = "SSD"
    }

    command <<<
        set -euxo pipefail

        mkdir -p ./tmp
        mkdir -p ./results

        export SPARK_LOCAL_DIRS="$PWD/tmp"

        DRIVER_MEM_GB=$((~{memory_gb} - 8))
        if [ "$DRIVER_MEM_GB" -lt 4 ]; then DRIVER_MEM_GB=4; fi

        export SPARK_DRIVER_MEMORY="${DRIVER_MEM_GB}g"
        export PYSPARK_SUBMIT_ARGS="--driver-memory ${DRIVER_MEM_GB}g --executor-memory ${DRIVER_MEM_GB}g pyspark-shell"
        export JAVA_OPTS="-Xms${DRIVER_MEM_GB}g -Xmx${DRIVER_MEM_GB}g"

        # Extract coverage.h5
        mkdir -p ./coverage_db
        tar -xzf ~{coverage_db_tar} -C ./coverage_db
        test -f ./coverage_db/coverage.h5

        # Extract input merged MT tar
        mkdir -p ./input_mt
        tar -xzf ~{in_mt_tar} -C ./input_mt
        IN_MT_DIR=$(find ./input_mt -maxdepth 2 -type d -name "*.mt" | head -n 1)
        if [ -z "$IN_MT_DIR" ]; then
            echo "ERROR: could not find .mt directory after extracting in_mt_tar" >&2
            find ./input_mt -maxdepth 3 -type d | head -100 >&2
            exit 1
        fi

        python3 /opt/mtSwirl/generate_mtdna_call_mt/Terra/finalize_mt_with_covdb.py \
            --in-mt "$IN_MT_DIR" \
            --coverage-h5-path ./coverage_db/coverage.h5 \
            --out-mt ./results/~{file_name}.mt \
            --temp-dir ./tmp \
            --minimum-homref-coverage ~{minimum_homref_coverage} \
            --homref-position-block-size ~{homref_position_block_size} \
            --artifact-prone-sites-path ~{artifact_prone_sites_path} \
            --artifact-prone-sites-reference ~{artifact_prone_sites_reference} \
            --n-final-partitions ~{n_final_partitions} \
            ~{if overwrite then "--overwrite" else ""}

        tar -czf results.tar.gz -C ./results "~{file_name}.mt"
    >>>

    output {
        File results_tar = "results.tar.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/aou-mitochondrial-combine-vcfs-covdb:dev-hailkeyfix"
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " " + disk_type
    }
}

task subset_data_table {
    input {
        File full_data_tsv
        File? sample_list_tsv

        # Runtime parameters
        Int memory_gb = 16
        Int cpu = 4 
        Int disk_gb = 100
    }

    String output_tsv = basename(select_first([sample_list_tsv, full_data_tsv]), ".tsv") + "_data.tsv"

    command <<<
    set -euxo pipefail
    python3 <<'EOF'
    import pandas as pd
    import sys

    # If sample_list_tsv is not defined, just copy the full TSV
    if "~{sample_list_tsv}" == "":
        df_main = pd.read_csv("~{full_data_tsv}", sep="\t", dtype=str)
        df_main.to_csv("~{output_tsv}", sep="\t", index=False)
        sys.exit(0)

    df_main = pd.read_csv("~{full_data_tsv}", sep="\t", dtype=str)
    df_samples = pd.read_csv("~{sample_list_tsv}", sep="\t", header=None, names=["sample_id"], dtype=str)

    # Check if TSV has header (Terra-style: entity:sample_id)
    first_col = df_main.columns[0]
    if first_col.startswith("entity:"):
        id_col = first_col
    else:
        sys.exit("ERROR: Unrecognized format for sample ID column in the full data TSV.")

    df_subset = df_main[df_main[id_col].isin(df_samples["sample_id"])]

    df_subset.to_csv("~{output_tsv}", sep="\t", index=False)
    EOF
    >>>

    output {
        File subset_tsv = "~{output_tsv}"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/python-numpy-pandas:1.0.0-2.2.3-1.25.2"
        memory: memory_gb + " GB" 
        cpu: cpu
        disks: "local-disk " + disk_gb + " HDD"
    }
}

task process_tsv_files {
    input {
        File coverage_tsv  # Path to genomics_metrics_Dec142023_1859_02_tz0000.tsv
        File ancestry_tsv  # Path to echo_v4_r2.ancestry_preds.tsv
        File dob_tsv       # Path to echo_DoB_data.tsv
        File wgs_median_coverage_tsv  # Path to wgs_median_coverage.tsv
        File input_tsv     # Input TSV file to process
        String output_tsv_name = "processed_data.tsv"  # Name of the output TSV file

        # Runtime parameters
        Int memory_gb = 16
        Int cpu = 4 
        Int disk_gb = 100
    }

    command <<<
        set -euxo pipefail

        python3 <<EOF

        import pandas as pd
        import numpy as np

        # Load the input TSV into a DataFrame
        df = pd.read_csv("~{input_tsv}", sep="\t")

        # Define only the required columns
        columns_needed = ["contamination", 
                          "coverage_metrics", 
                          "final_base_level_coverage_metrics", 
                          "major_haplogroup",
                          "mean_coverage", 
                          "median_coverage", 
                          "mtdna_consensus_overlaps",
                          "final_vcf", 
                          "research_id"] 
        
        # Keep only necessary columns
        #filtered_df = df[columns_needed]
     
        # We're actually just going to keep all the columns and see if that fixes a type issue that is arising downstream
        filtered_df = df

        # Add and rename columns
        filtered_df["s"] = filtered_df["research_id"]
        filtered_df = filtered_df.rename(columns={
            "research_id": "entity:participant_id",
            "final_base_level_coverage_metrics": "coverage",
            "mean_coverage": "mt_mean_coverage",
            "median_coverage": "mt_median_coverage"
        })

        # Load additional TSV files
        coveragetsv_df = pd.read_csv("~{coverage_tsv}", sep="\t")
        ancestrytsv_df = pd.read_csv("~{ancestry_tsv}", sep="\t")
        dobtsv_df = pd.read_csv("~{dob_tsv}", sep="\t")
        medcoveragetsv_df = pd.read_csv("~{wgs_median_coverage_tsv}", sep="\t")

        # Merge with filtered_df on 'research_id' and 's'
        filtered_df = filtered_df.merge(
            coveragetsv_df[['research_id', 'mean_coverage', 'biosample_collection_date', 'verify_bam_id2_contamination']],
            left_on='s', right_on='research_id', how='left'
        )
        filtered_df.drop(columns=['research_id'], inplace=True)
        filtered_df = filtered_df.merge(
            ancestrytsv_df[['research_id', 'ancestry_pred']],
            left_on='s', right_on='research_id', how='left'
        )
        filtered_df.drop(columns=['research_id'], inplace=True)
        filtered_df = filtered_df.merge(
            dobtsv_df[['research_id', 'date_of_birth']],
            left_on='s', right_on='research_id', how='left'
        )
        filtered_df.drop(columns=['research_id'], inplace=True)
        filtered_df = filtered_df.merge(
            medcoveragetsv_df[['research_id', 'median_coverage']],
            left_on='s', right_on='research_id', how='left'
        )
        filtered_df.drop(columns=['research_id'], inplace=True)

        # Add a check to ensure that our filtered_df has the same number of samples (shape[0]) as the original df
        if filtered_df.shape[0] != df.shape[0]:
            raise ValueError("Filtered DataFrame does not have the same number of samples as the original.")

        # Calculate age
        filtered_df['date_of_birth'] = pd.to_datetime(filtered_df['date_of_birth'])
        filtered_df['biosample_collection_date'] = pd.to_datetime(filtered_df['biosample_collection_date'])
        filtered_df['age'] = pd.to_numeric(
            np.floor((filtered_df['biosample_collection_date'] - filtered_df['date_of_birth']).dt.days / 365)
        )

        # Age must be an int and must be present
        filtered_df['age'] = filtered_df['age'].astype(int)
        if filtered_df['age'].isna().any():
            raise ValueError("Unexpected missing ages detected.")

        # Rename columns for compatibility
        filtered_df.rename(columns={"mean_coverage": "wgs_mean_coverage"}, inplace=True)
        filtered_df.rename(columns={"median_coverage": "wgs_median_coverage"}, inplace=True)
        filtered_df.rename(columns={"ancestry_pred": "pop"}, inplace=True)
        filtered_df.rename(columns={"verify_bam_id2_contamination": "freemix_percentage"}, inplace=True)

        # Filter rows with valid coverage metrics
        filtered_df = filtered_df[
            (filtered_df['coverage_metrics'].notna()) & (filtered_df['coverage_metrics'] != '')
        ]

        # Save the processed DataFrame to a TSV file
        filtered_df.to_csv("~{output_tsv_name}", sep="\t", index=False)
        EOF
    >>>

    output {
        File processed_tsv = "~{output_tsv_name}"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/python-numpy-pandas:1.0.0-2.2.3-1.25.2"
        memory: memory_gb + " GB" 
        cpu: cpu
        disks: "local-disk " + disk_gb + " HDD"
    }
}

task annotate_coverage {
    input {
        File input_tsv        # Input TSV file
        # v2 coverage DB builder parameters (kept in this task for drop-in workflow wiring)
        Int batch_size = 256
        Int position_block_size = 4096

        Boolean skip_summary = false

        # Runtime parameters
        Int memory_gb = 256 #1000
        Int cpu = 16 #64
        Int disk_gb = 2000 #2000
        String disk_type = "SSD" #"SSD"
    }

    command <<<
        set -euxo pipefail

        WORK_DIR=$(pwd)
        # Build HDF5 coverage DB + summary TSV.
        # Note: this image contains only the coverage_db builder + deps (no Hail/Spark).
        python -m generate_mtdna_call_mt.coverage_db.build_coverage_db \
            --input-tsv ~{input_tsv} \
            --out-h5 coverage.h5 \
            ~{if skip_summary then "--skip-summary" else "--out-summary-tsv coverage_summary.tsv"} \
            --batch-size ~{batch_size} \
            --position-block-size ~{position_block_size}

        ls -lh coverage.h5
        if [ "~{skip_summary}" = "true" ]; then
            echo "skip_summary enabled: not generating coverage_summary.tsv"
            # Bundle outputs as a single artifact for consistent WDL handling.
            tar -czf $WORK_DIR/coverage_db.tar.gz coverage.h5
        else
            ls -lh coverage_summary.tsv
            # Bundle outputs as a single artifact for consistent WDL handling.
            tar -czf $WORK_DIR/coverage_db.tar.gz coverage.h5 coverage_summary.tsv
        fi
    >>>

    output {
        File output_ht = "coverage_db.tar.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/mt-coverage-db:dev"
        memory: memory_gb + " GB" 
        cpu: cpu
        disks: "local-disk " + disk_gb + " " + disk_type 
    }
}

task combine_vcfs {
    input {
        File input_tsv        # Input TSV file
        File coverage_mt_tar    # Tar.gzipped directory of the Hail table
        String a_ref = "GRCh38"           # Reference genome (e.g., GRCh38)
        Boolean overwrite = false  # Overwrite existing files (default: false)
        String vcf_col_name = "final_vcf"     # Column name for VCFs
        String artifact_prone_sites_path = "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed"  # Path to artifact-prone sites BED file
        String file_name        # Output file name

        # Runtime parameters
        Int memory_gb = 64 #1000
        Int cpu = 16 #64
        Int disk_gb = 500 #2000
        String disk_type = "HDD" #"SSD"
    }

    command <<<
        set -euxo pipefail

        mkdir -p ./tmp
        mkdir -p ./results

        WORK_DIR=$(pwd)

        # Set the spark config via environment variables
        export SPARK_LOCAL_DIRS="$PWD/tmp"

        DRIVER_MEM_GB=$((~{memory_gb} - 4))
        if [ "$DRIVER_MEM_GB" -lt 2 ]; then DRIVER_MEM_GB=4; fi

        export SPARK_DRIVER_MEMORY="${DRIVER_MEM_GB}g"
        export PYSPARK_SUBMIT_ARGS="--driver-memory ${DRIVER_MEM_GB}g --executor-memory ${DRIVER_MEM_GB}g pyspark-shell"
        export JAVA_OPTS="-Xms${DRIVER_MEM_GB}g -Xmx${DRIVER_MEM_GB}g"  

        # Unzip the tar.gz file containing the Hail table
        mkdir -p ./unzipped_coverage.ht
        tar -xzf ~{coverage_mt_tar} -C /cromwell_root/unzipped_coverage.ht
        ls -lh /cromwell_root/unzipped_coverage.ht
        cp -r /cromwell_root/unzipped_coverage.ht/merged_coverage_tsvs.mt /cromwell_root/unzipped_coverage.mt

        # Verify the extracted directory TODO
        if [ ! -d "/cromwell_root/unzipped_coverage.mt" ]; then
            echo "Error: Directory '/cromwell_root/unzipped_coverage.mt' does not exist after extraction."
            exit 1
        fi


        # Run the combine_vcfs.py script
        python3 /opt/mtSwirl/generate_mtdna_call_mt/Terra/combine_vcfs.py \
        --input-tsv ~{input_tsv} \
        -c /cromwell_root/unzipped_coverage.mt \
        -a-ref ~{a_ref} \
    ~{if overwrite then "--overwrite" else ""} \
        --vcf-col-name ~{vcf_col_name} \
        --output-bucket ./results \
        --temp-dir ./tmp \
        --artifact-prone-sites-path ~{artifact_prone_sites_path} \
        --file-name ~{file_name} \
        --include-extra-v2-fields


        # Tar zip the results directory
        tar -czf $WORK_DIR/results.tar.gz ./results/combined_vcf.vcf.gz.mt
    >>>

    output {
        File results_tar = "results.tar.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/aou-mitochondrial-annotate-coverage:1.0.0"
        memory: memory_gb + " GB" 
        cpu: cpu
        disks: "local-disk " + disk_gb + " " + disk_type 
    }
}

task combine_vcfs_and_homref_from_covdb {
    input {
        File input_tsv              # Input TSV file; contains gs:// VCFs in vcf_col_name
        File coverage_db_tar        # Tar.gz produced by annotate_coverage containing coverage.h5 (+ optional summary)
        String vcf_col_name = "final_vcf"
        String artifact_prone_sites_path = "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed"
        String artifact_prone_sites_reference = "default"
        String file_name            # Output base name (writes <file_name>.mt)
        Int minimum_homref_coverage = 100
        Int chunk_size = 100
        Int homref_position_block_size = 1024
        Int n_final_partitions = 1000
        Int split_merging = 1
        Boolean overwrite = false
        Boolean include_extra_v2_fields = true

        # Runtime parameters
        Int memory_gb = 256
        Int cpu = 32
        Int disk_gb = 2000
        String disk_type = "SSD"
    }

    command <<<
        set -euxo pipefail

        mkdir -p ./tmp
        mkdir -p ./results

        WORK_DIR=$(pwd)

        # Spark/Hail scratch
        export SPARK_LOCAL_DIRS="$PWD/tmp"

        DRIVER_MEM_GB=$((~{memory_gb} - 8))
        if [ "$DRIVER_MEM_GB" -lt 4 ]; then DRIVER_MEM_GB=4; fi

        export SPARK_DRIVER_MEMORY="${DRIVER_MEM_GB}g"
        export PYSPARK_SUBMIT_ARGS="--driver-memory ${DRIVER_MEM_GB}g --executor-memory ${DRIVER_MEM_GB}g pyspark-shell"
        export JAVA_OPTS="-Xms${DRIVER_MEM_GB}g -Xmx${DRIVER_MEM_GB}g"

        # Extract coverage.h5 from tarball
        mkdir -p ./coverage_db
        tar -xzf ~{coverage_db_tar} -C ./coverage_db
        ls -lh ./coverage_db
        test -f ./coverage_db/coverage.h5

        python3 /opt/mtSwirl/generate_mtdna_call_mt/Terra/combine_vcfs_and_homref_from_covdb.py \
            --input-tsv ~{input_tsv} \
            --coverage-h5-path ./coverage_db/coverage.h5 \
            --vcf-col-name ~{vcf_col_name} \
            --artifact-prone-sites-path ~{artifact_prone_sites_path} \
            --artifact-prone-sites-reference ~{artifact_prone_sites_reference} \
            --output-bucket ./results \
            --temp-dir ./tmp \
            --file-name ~{file_name} \
            --minimum-homref-coverage ~{minimum_homref_coverage} \
            --chunk-size ~{chunk_size} \
            --homref-position-block-size ~{homref_position_block_size} \
            --n-final-partitions ~{n_final_partitions} \
            --split-merging ~{split_merging} \
            ~{if overwrite then "--overwrite" else ""} \
            ~{if include_extra_v2_fields then "--include-extra-v2-fields" else ""}

        # Tar the MT directory.
        # The script writes: ./results/<file_name>.mt
        tar -czf $WORK_DIR/results.tar.gz -C ./results "~{file_name}.mt"
    >>>

    output {
        File results_tar = "results.tar.gz"
    }

    runtime {
        # NOTE: This must be a Hail-capable image with mtSwirl code baked in at /opt/mtSwirl.
        docker: "us.gcr.io/broad-gotc-prod/aou-mitochondrial-combine-vcfs-covdb:dev-hailkeyfix"
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " " + disk_type
    }
}

task add_annotations {
    input {
        File coverage_mt_tar    # Tar.gzipped directory of the Hail table
        Boolean keep_all_samples = false  # Keep all samples (default: false)
        File coverage_tsv     # Path to the coverage input TSV file
        File vcf_mt             # Path to the MatrixTable
        String output_name      # directory output name
        
        # Runtime parameters
        Int memory_gb = 64 #1000
        Int cpu = 16 #64
        Int disk_gb = 500 #2000
        String disk_type = "HDD" #"SSD"
    }

     command <<<
        set -euxo pipefail

        WORK_DIR=$(pwd)

        # Set the spark config via environment variables
        export SPARK_LOCAL_DIRS="$PWD/tmp"

        DRIVER_MEM_GB=$((~{memory_gb} - 4))
        if [ "$DRIVER_MEM_GB" -lt 2 ]; then DRIVER_MEM_GB=4; fi

        export SPARK_DRIVER_MEMORY="${DRIVER_MEM_GB}g"
        export PYSPARK_SUBMIT_ARGS="--driver-memory ${DRIVER_MEM_GB}g --executor-memory ${DRIVER_MEM_GB}g pyspark-shell"
        export JAVA_OPTS="-Xms${DRIVER_MEM_GB}g -Xmx${DRIVER_MEM_GB}g"  

        # Unzip coverage MatrixTable tarball
        mkdir -p ./unzipped_coverage.mt
        tar -xzf ~{coverage_mt_tar} -C ./unzipped_coverage.mt
        ls -lh ./unzipped_coverage.mt/merged_coverage_tsvs.mt

        # Unzip VCF MatrixTable tarball
        mkdir -p ./unzipped_vcf.mt
        tar -xzf ~{vcf_mt} -C ./unzipped_vcf.mt
        ls -lh ./unzipped_vcf.mt/results/combined_vcf.vcf.gz.mt

        # Verify extraction
        if [ ! -d "./unzipped_coverage.mt" ]; then
            echo "Error: Directory './unzipped_coverage.mt' does not exist after extraction."
            exit 1
        fi

        # Run the add_annotations.py script baked inside mtSwirl clone
        python3 /opt/mtSwirl/generate_mtdna_call_mt/add_annotations.py \
            --sample-stats=~{coverage_tsv} \
            ~{if keep_all_samples then "--keep-all-samples" else ""} \
            --fully-skip-vep \
            --band-aid-dbsnp-path-fix \
            --min-het-threshold 0.05 \
            -v ./~{output_name}/vep \
            -a ~{coverage_tsv} \
            -m unzipped_vcf.mt/results/combined_vcf.vcf.gz.mt \
            -d ./~{output_name}

        # Compress the annotated output directory
        tar -czf $WORK_DIR/annotated_output.tar.gz ~{output_name}
    >>>

    output {
        File annotated_output_tar = "annotated_output.tar.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/aou-mitochondrial-annotate-coverage:1.0.0"
        memory: memory_gb + " GB" 
        cpu: cpu
        disks: "local-disk " + disk_gb + " " + disk_type 
    }
}