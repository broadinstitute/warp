version 1.0

workflow mt_coverage_merge {

    input {
        # Side inputs with fields that are not incoporated in the samples table
        String coverage_tsv
        String ancestry_tsv
        String dob_tsv


        File full_data_tsv
        File? sample_list_tsv

    }

    String pipeline_version = "1.0.0"

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
            input_tsv = input_table
    }

    call annotate_coverage {
        input:
            input_tsv = process_tsv_files.processed_tsv  # Input TSV file path
    }

    call combine_vcfs {
        input:
            input_tsv = process_tsv_files.processed_tsv,  # Input TSV file path
            coverage_mt_tar = annotate_coverage.output_ht,  # Tar.gzipped directory of the Hail table
            artifact_prone_sites_path = "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed",  # Path to artifact-prone sites BED file
            file_name = "combined_vcf.vcf.gz"  # Output file name
    }

    call add_annotations as annotated {
        input:
            coverage_mt_tar = annotate_coverage.output_ht,  # Tar.gzipped directory of the Hail table
            coverage_tsv = process_tsv_files.processed_tsv,  # Path to the coverage input TSV file
            vcf_mt = combine_vcfs.results_tar,  # Path to the MatrixTable
            keep_all_samples = true,
            output_name = "annotated"
    }

    call add_annotations as filt_annotated {
        input:
            coverage_mt_tar = annotate_coverage.output_ht,  # Tar.gzipped directory of the Hail table
            coverage_tsv = process_tsv_files.processed_tsv,  # Path to the coverage input TSV file
            vcf_mt = combine_vcfs.results_tar,  # Path to the MatrixTable
            keep_all_samples = false,
            output_name = "filt_annotated"
    }

    output {
        File processed_tsv = process_tsv_files.processed_tsv
        File output_coverage_ht = annotate_coverage.output_ht
        File combined_vcf = combine_vcfs.results_tar
        File annotated_output_tar = annotated.annotated_output_tar
        File filt_annotated_output_tar = filt_annotated.annotated_output_tar
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
                          "mt_final_vcf", 
                          "research_id",
                          "sample"] 
        
        # Keep only necessary columns
        filtered_df = df[columns_needed]

        # Add and rename columns
        filtered_df["s"] = filtered_df["research_id"]
        filtered_df = filtered_df.rename(columns={
            "research_id": "entity:participant_id",
            "final_base_level_coverage_metrics": "coverage",
            "mean_coverage": "mt_mean_coverage",
            "median_coverage": "mt_median_coverage",
            "mt_final_vcf": "final_vcf"
        })

        # Load additional TSV files
        coveragetsv_df = pd.read_csv("~{coverage_tsv}", sep="\t")
        ancestrytsv_df = pd.read_csv("~{ancestry_tsv}", sep="\t")
        dobtsv_df = pd.read_csv("~{dob_tsv}", sep="\t")

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

        # Add a check to ensure that our filtered_df has the same number of samples (shape[0]) as the original df
        if filtered_df.shape[0] != df.shape[0]:
            raise ValueError("Filtered DataFrame does not have the same number of samples as the original.")

        # Calculate age
        filtered_df['date_of_birth'] = pd.to_datetime(filtered_df['date_of_birth'])
        filtered_df['biosample_collection_date'] = pd.to_datetime(filtered_df['biosample_collection_date'])
        filtered_df['age'] = pd.to_numeric(
            np.floor((filtered_df['biosample_collection_date'] - filtered_df['date_of_birth']).dt.days / 365)
        )
        # commenting this out, since we want to test what happens with missing ages
        #filtered_df['age'] = filtered_df['age'].fillna(39).astype(int)

        # we will rename the verify_bam_id2_contamination to freemix_percentage instead of setting to 0 
        #filtered_df["freemix_percentage"] = 0 

        # Rename columns for compatibility
        filtered_df.rename(columns={"mean_coverage": "wgs_mean_coverage"}, inplace=True)
        filtered_df.rename(columns={"ancestry_pred": "pop"}, inplace=True)
        filtered_df.rename(columns={"verify_bam_id2_contamination": "freemix_percentage"}, inplace=True)


        # Temporary workaround - we are working on getting the real values from the Dragen metrics
        filtered_df["wgs_median_coverage"] = filtered_df["wgs_mean_coverage"]

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
        String input_tsv        # Input TSV file path (optional)
        Int? chunk_size = 100    # Chunk size for combining VCFs (default: 100)
        Boolean overwrite = false  # Overwrite existing files (default: false)
        Boolean keep_targets = false  # Add annotation for target (default: false)
        Boolean hail_only = false  # Skip generating flat files (default: false)
        Int? split_merging = 10  # Number of jobs for splitting merging (default: 1)

        # Runtime parameters
        Int memory_gb = 32 #1000
        Int cpu = 4 #64
        Int disk_gb = 500 #2000
        String disk_type = "HDD" #"SSD"
        String cpu_platform = "Intel Ice Lake"
    }

    command <<<
        set -euxo pipefail

        mkdir -p ./tmp
        mkdir -p ./results.ht

        WORK_DIR=$(pwd)


        # Run the annotate_coverage.py script
        python3 /opt/mtSwirl/generate_mtdna_call_mt/Terra/annotate_coverage.py \
        ~{if overwrite then "--overwrite" else ""} \
        ~{if keep_targets then "--keep-targets" else ""} \
        --input-tsv ~{input_tsv} \
        --output-ht "./merged_coverage_tsvs.ht" \
        --temp-dir "./tmp/" \
        --chunk-size ~{chunk_size} \
        ~{if hail_only then "--hail-only" else ""} \
        --split-merging ~{split_merging}


        ## note that both the ht and the mt are outputted by the tool
        ls -lh ./merged_coverage_tsvs.ht
        ls -lht ./merged_coverage_tsvs.mt

        # Archive the MatrixTable
        tar -czf $WORK_DIR/coverages_tsv.mt.tar.gz ./merged_coverage_tsvs.mt*
    >>>

    output {
        File output_ht = "coverages_tsv.mt.tar.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/aou-mitochondrial-annotate-coverage:1.0.0"
        memory: memory_gb + " GB" 
        cpu: cpu
        disks: "local-disk " + disk_gb + " " + disk_type 
        cpu_platform: cpu_platform
    }
}

task combine_vcfs {
    input {
        String input_tsv        # Input TSV file path
        File coverage_mt_tar    # Tar.gzipped directory of the Hail table
        String a_ref = "GRCh38"           # Reference genome (e.g., GRCh38)
        Boolean overwrite = false  # Overwrite existing files (default: false)
        String vcf_col_name = "final_vcf"     # Column name for VCFs
        String output_bucket = "./results"    # Output bucket path
        String artifact_prone_sites_path = "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed"  # Path to artifact-prone sites BED file
        String file_name        # Output file name

        # Runtime parameters
        Int memory_gb = 32 #1000
        Int cpu = 4 #64
        Int disk_gb = 500 #2000
        String disk_type = "HDD" #"SSD"
        String cpu_platform = "Intel Ice Lake"
    }

    command <<<
        set -euxo pipefail

        mkdir -p ./tmp
        mkdir -p ./results

        WORK_DIR=$(pwd)

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
        --overwrite \
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
        cpu_platform: cpu_platform
    }
}

task add_annotations {
    input {
        File coverage_mt_tar    # Tar.gzipped directory of the Hail table
        Boolean keep_all_samples = false  # Keep all samples (default: false)
        String coverage_tsv     # Path to the coverage input TSV file
        File vcf_mt             # Path to the MatrixTable
        String output_name      # directory output name
        
        # Runtime parameters
        Int memory_gb = 32 #1000
        Int cpu = 4 #64
        Int disk_gb = 500 #2000
        String disk_type = "HDD" #"SSD"
        String cpu_platform = "Intel Ice Lake"
    }

     command <<<
        set -euxo pipefail

        WORK_DIR=$(pwd)

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
        cpu_platform: cpu_platform
    }
}