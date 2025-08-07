version 1.0

workflow mt_coverage_merge {

    input {
        # Side inputs with fields that are not incoporated in the samples table
        String coverage_tsv
        String ancestry_tsv
        String dob_tsv

        String hail_docker_img = "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-hail:1.1.12"

        ## Copy of all of the fields in the sample table (this can be reduced)
        Array[String] n_liftover_r2_spanning_complex
        Array[String] success_liftover_variants
        Array[String] stats_outputs
        Array[String] self_mt_aligned_bai
        Array[String] n_liftover_r2_repaired_success
        Array[String] r1_split_vcf_index
        Array[String] cram_intact
        Array[String] self_base_level_coverage_metrics
        Array[String] self_reference_fasta
        Array[String] research_id
        Array[String] subset_bam
        Array[String] self_ref_split_vcf_index
        Array[String] r1_nuc_vcf_index
        Array[String] n_liftover_r2_failed_het_dele_span_insertion_boundary
        Array[String] self_control_region_shifted
        Array[String] mt_final_rejected_vcf
        Array[String] contamination_metrics
        Array[String] selfSM
        Array[String] nuc_variants_pass
        Array[String] n_liftover_r2_spanningfixrhs_sharedlhs
        Array[String] reblocked_gvcf_index
        Array[String] n_liftover_r2_ref_insertion_new_haplo
        Array[String] r1_nuc_vcf_unfiltered
        Array[String] mt_median_coverage
        Array[String] n_liftover_r2_failed_new_dupes_leftshift
        Array[String] r1_nuc_vcf
        Array[String] input_vcf_for_haplochecker
        Array[String] subset_bai
        Array[String] n_liftover_r2_injected_from_success
        Array[String] reblocked_gvcf
        Array[String] crai_path
        Array[String] numt_base_level_coverage
        Array[String] self_ref_vcf
        Array[String] self_mt_aligned_bam
        Array[String] coverage_metrics
        Array[String] fixed_liftover_variants
        Array[String] reference_to_self_ref_chain
        Array[String] r1_split_vcf
        Array[String] self_ref_vcf_index
        Array[String] self_ref_split_vcf
        Array[String] contamination
        Array[String] self_non_control_region
        Array[String] failed_liftover_variants
        Array[String] mt_final_vcf
        Array[String] n_reads_unpaired_dropped
        Array[String] liftover_fix_pipeline_log
        Array[String] major_haplogroup
        Array[String] passes_qc
        Array[String] n_liftover_r2_het_ins_sharing_lhs_hom_dele
        Array[String] n_liftover_r2_left_shift
        Array[String] final_base_level_coverage_metrics
        Array[String] sex_at_birth
        Array[String] mt_mean_coverage
        Array[String] r1_vcf_index
        Array[String] n_liftover_r2_spanningfixlhs_upstream
        Array[String] nuc_consensus_overlaps
        Array[String] mtdna_consensus_overlaps
        Array[String] nuc_variants_dropped
        Array[String] r1_vcf
        Array[String] duplicate_metrics
        Array[String] theoretical_sensitivity_metrics
        Array[String] cram_path
    }

    call generate_tsv_file {
        input:
            n_liftover_r2_spanning_complex = n_liftover_r2_spanning_complex,
            success_liftover_variants = success_liftover_variants,
            stats_outputs = stats_outputs,
            self_mt_aligned_bai = self_mt_aligned_bai,
            n_liftover_r2_repaired_success = n_liftover_r2_repaired_success,
            r1_split_vcf_index = r1_split_vcf_index,
            cram_intact = cram_intact,
            self_base_level_coverage_metrics = self_base_level_coverage_metrics,
            self_reference_fasta = self_reference_fasta,
            research_id = research_id,
            subset_bam = subset_bam,
            self_ref_split_vcf_index = self_ref_split_vcf_index,
            r1_nuc_vcf_index = r1_nuc_vcf_index,
            n_liftover_r2_failed_het_dele_span_insertion_boundary = n_liftover_r2_failed_het_dele_span_insertion_boundary,
            self_control_region_shifted = self_control_region_shifted,
            mt_final_rejected_vcf = mt_final_rejected_vcf,
            contamination_metrics = contamination_metrics,
            selfSM = selfSM,
            nuc_variants_pass = nuc_variants_pass,
            n_liftover_r2_spanningfixrhs_sharedlhs = n_liftover_r2_spanningfixrhs_sharedlhs,
            reblocked_gvcf_index = reblocked_gvcf_index,
            n_liftover_r2_ref_insertion_new_haplo = n_liftover_r2_ref_insertion_new_haplo,
            r1_nuc_vcf_unfiltered = r1_nuc_vcf_unfiltered,
            mt_median_coverage = mt_median_coverage,
            n_liftover_r2_failed_new_dupes_leftshift = n_liftover_r2_failed_new_dupes_leftshift,
            r1_nuc_vcf = r1_nuc_vcf,
            input_vcf_for_haplochecker = input_vcf_for_haplochecker,
            subset_bai = subset_bai,
            n_liftover_r2_injected_from_success = n_liftover_r2_injected_from_success,
            reblocked_gvcf = reblocked_gvcf,
            crai_path = crai_path,
            numt_base_level_coverage = numt_base_level_coverage,
            self_ref_vcf = self_ref_vcf,
            self_mt_aligned_bam = self_mt_aligned_bam,
            coverage_metrics = coverage_metrics,
            fixed_liftover_variants = fixed_liftover_variants,
            reference_to_self_ref_chain = reference_to_self_ref_chain,
            r1_split_vcf = r1_split_vcf,
            self_ref_vcf_index = self_ref_vcf_index,
            self_ref_split_vcf = self_ref_split_vcf,
            contamination = contamination,
            self_non_control_region = self_non_control_region,
            failed_liftover_variants = failed_liftover_variants,
            mt_final_vcf = mt_final_vcf,
            n_reads_unpaired_dropped = n_reads_unpaired_dropped,
            liftover_fix_pipeline_log = liftover_fix_pipeline_log,
            major_haplogroup = major_haplogroup,
            passes_qc = passes_qc,
            n_liftover_r2_het_ins_sharing_lhs_hom_dele = n_liftover_r2_het_ins_sharing_lhs_hom_dele,
            n_liftover_r2_left_shift = n_liftover_r2_left_shift,
            final_base_level_coverage_metrics = final_base_level_coverage_metrics,
            sex_at_birth = sex_at_birth,
            mt_mean_coverage = mt_mean_coverage,
            r1_vcf_index = r1_vcf_index,
            n_liftover_r2_spanningfixlhs_upstream = n_liftover_r2_spanningfixlhs_upstream,
            nuc_consensus_overlaps = nuc_consensus_overlaps,
            mtdna_consensus_overlaps = mtdna_consensus_overlaps,
            nuc_variants_dropped = nuc_variants_dropped,
            r1_vcf = r1_vcf,
            duplicate_metrics = duplicate_metrics,
            theoretical_sensitivity_metrics = theoretical_sensitivity_metrics,
            cram_path = cram_path
    }

    call process_tsv_files {
        input:
            coverage_tsv = coverage_tsv,
            ancestry_tsv = ancestry_tsv,
            dob_tsv = dob_tsv,
            input_tsv = generate_tsv_file.summary_tsv,
    }

    call annotate_coverage {
        input:
            input_tsv = process_tsv_files.processed_tsv,  # Input TSV file path
            hail_docker = hail_docker_img  # Docker image with Hail and dependencies
    }

    output {
        File processed_tsv = process_tsv_files.processed_tsv
        #        Directory output_ht = annotate_coverage.output_ht
        File output_ht = annotate_coverage.output_ht
    }
}

task generate_tsv_file {

    input {
        ## Copy inputs from the matrix table... this is overkill
        Array[String] n_liftover_r2_spanning_complex
        Array[String] success_liftover_variants
        Array[String] stats_outputs
        Array[String] self_mt_aligned_bai
        Array[String] n_liftover_r2_repaired_success
        Array[String] r1_split_vcf_index
        Array[String] cram_intact
        Array[String] self_base_level_coverage_metrics
        Array[String] self_reference_fasta
        Array[String] research_id
        Array[String] subset_bam
        Array[String] self_ref_split_vcf_index
        Array[String] r1_nuc_vcf_index
        Array[String] n_liftover_r2_failed_het_dele_span_insertion_boundary
        Array[String] self_control_region_shifted
        Array[String] mt_final_rejected_vcf
        Array[String] contamination_metrics
        Array[String] selfSM
        Array[String] nuc_variants_pass
        Array[String] n_liftover_r2_spanningfixrhs_sharedlhs
        Array[String] reblocked_gvcf_index
        Array[String] n_liftover_r2_ref_insertion_new_haplo
        Array[String] r1_nuc_vcf_unfiltered
        Array[String] mt_median_coverage
        Array[String] n_liftover_r2_failed_new_dupes_leftshift
        Array[String] r1_nuc_vcf
        Array[String] input_vcf_for_haplochecker
        Array[String] subset_bai
        Array[String] n_liftover_r2_injected_from_success
        Array[String] reblocked_gvcf
        Array[String] crai_path
        Array[String] numt_base_level_coverage
        Array[String] self_ref_vcf
        Array[String] self_mt_aligned_bam
        Array[String] coverage_metrics
        Array[String] fixed_liftover_variants
        Array[String] reference_to_self_ref_chain
        Array[String] r1_split_vcf
        Array[String] self_ref_vcf_index
        Array[String] self_ref_split_vcf
        Array[String] contamination
        Array[String] self_non_control_region
        Array[String] failed_liftover_variants
        Array[String] mt_final_vcf
        Array[String] n_reads_unpaired_dropped
        Array[String] liftover_fix_pipeline_log
        Array[String] major_haplogroup
        Array[String] passes_qc
        Array[String] n_liftover_r2_het_ins_sharing_lhs_hom_dele
        Array[String] n_liftover_r2_left_shift
        Array[String] final_base_level_coverage_metrics
        Array[String] sex_at_birth
        Array[String] mt_mean_coverage
        Array[String] r1_vcf_index
        Array[String] n_liftover_r2_spanningfixlhs_upstream
        Array[String] nuc_consensus_overlaps
        Array[String] mtdna_consensus_overlaps
        Array[String] nuc_variants_dropped
        Array[String] r1_vcf
        Array[String] duplicate_metrics
        Array[String] theoretical_sensitivity_metrics
        Array[String] cram_path
    }

    command <<<
        # Write all input arrays to a temporary file
        echo -e "n_liftover_r2_spanning_complex\t~{sep="\t" n_liftover_r2_spanning_complex}" > input_data.tsv
        echo -e "success_liftover_variants\t~{sep="\t" success_liftover_variants}" >> input_data.tsv
        echo -e "stats_outputs\t~{sep="\t" stats_outputs}" >> input_data.tsv
        echo -e "self_mt_aligned_bai\t~{sep="\t" self_mt_aligned_bai}" >> input_data.tsv
        echo -e "n_liftover_r2_repaired_success\t~{sep="\t" n_liftover_r2_repaired_success}" >> input_data.tsv
        echo -e "r1_split_vcf_index\t~{sep="\t" r1_split_vcf_index}" >> input_data.tsv
        echo -e "cram_intact\t~{sep="\t" cram_intact}" >> input_data.tsv
        echo -e "self_base_level_coverage_metrics\t~{sep="\t" self_base_level_coverage_metrics}" >> input_data.tsv
        echo -e "self_reference_fasta\t~{sep="\t" self_reference_fasta}" >> input_data.tsv
        echo -e "research_id\t~{sep="\t" research_id}" >> input_data.tsv
        echo -e "subset_bam\t~{sep="\t" subset_bam}" >> input_data.tsv
        echo -e "self_ref_split_vcf_index\t~{sep="\t" self_ref_split_vcf_index}" >> input_data.tsv
        echo -e "r1_nuc_vcf_index\t~{sep="\t" r1_nuc_vcf_index}" >> input_data.tsv
        echo -e "n_liftover_r2_failed_het_dele_span_insertion_boundary\t~{sep="\t" n_liftover_r2_failed_het_dele_span_insertion_boundary}" >> input_data.tsv
        echo -e "self_control_region_shifted\t~{sep="\t" self_control_region_shifted}" >> input_data.tsv
        echo -e "mt_final_rejected_vcf\t~{sep="\t" mt_final_rejected_vcf}" >> input_data.tsv
        echo -e "contamination_metrics\t~{sep="\t" contamination_metrics}" >> input_data.tsv
        echo -e "selfSM\t~{sep="\t" selfSM}" >> input_data.tsv
        echo -e "nuc_variants_pass\t~{sep="\t" nuc_variants_pass}" >> input_data.tsv
        echo -e "n_liftover_r2_spanningfixrhs_sharedlhs\t~{sep="\t" n_liftover_r2_spanningfixrhs_sharedlhs}" >> input_data.tsv
        echo -e "reblocked_gvcf_index\t~{sep="\t" reblocked_gvcf_index}" >> input_data.tsv
        echo -e "n_liftover_r2_ref_insertion_new_haplo\t~{sep="\t" n_liftover_r2_ref_insertion_new_haplo}" >> input_data.tsv
        echo -e "r1_nuc_vcf_unfiltered\t~{sep="\t" r1_nuc_vcf_unfiltered}" >> input_data.tsv
        echo -e "mt_median_coverage\t~{sep="\t" mt_median_coverage}" >> input_data.tsv
        echo -e "n_liftover_r2_failed_new_dupes_leftshift\t~{sep="\t" n_liftover_r2_failed_new_dupes_leftshift}" >> input_data.tsv
        echo -e "r1_nuc_vcf\t~{sep="\t" r1_nuc_vcf}" >> input_data.tsv
        echo -e "input_vcf_for_haplochecker\t~{sep="\t" input_vcf_for_haplochecker}" >> input_data.tsv
        echo -e "subset_bai\t~{sep="\t" subset_bai}" >> input_data.tsv
        echo -e "n_liftover_r2_injected_from_success\t~{sep="\t" n_liftover_r2_injected_from_success}" >> input_data.tsv
        echo -e "reblocked_gvcf\t~{sep="\t" reblocked_gvcf}" >> input_data.tsv
        echo -e "crai_path\t~{sep="\t" crai_path}" >> input_data.tsv
        echo -e "numt_base_level_coverage\t~{sep="\t" numt_base_level_coverage}" >> input_data.tsv
        echo -e "self_ref_vcf\t~{sep="\t" self_ref_vcf}" >> input_data.tsv
        echo -e "self_mt_aligned_bam\t~{sep="\t" self_mt_aligned_bam}" >> input_data.tsv
        echo -e "coverage_metrics\t~{sep="\t" coverage_metrics}" >> input_data.tsv
        echo -e "fixed_liftover_variants\t~{sep="\t" fixed_liftover_variants}" >> input_data.tsv
        echo -e "reference_to_self_ref_chain\t~{sep="\t" reference_to_self_ref_chain}" >> input_data.tsv
        echo -e "r1_split_vcf\t~{sep="\t" r1_split_vcf}" >> input_data.tsv
        echo -e "self_ref_vcf_index\t~{sep="\t" self_ref_vcf_index}" >> input_data.tsv
        echo -e "self_ref_split_vcf\t~{sep="\t" self_ref_split_vcf}" >> input_data.tsv
        echo -e "contamination\t~{sep="\t" contamination}" >> input_data.tsv
        echo -e "self_non_control_region\t~{sep="\t" self_non_control_region}" >> input_data.tsv
        echo -e "failed_liftover_variants\t~{sep="\t" failed_liftover_variants}" >> input_data.tsv
        echo -e "mt_final_vcf\t~{sep="\t" mt_final_vcf}" >> input_data.tsv
        echo -e "n_reads_unpaired_dropped\t~{sep="\t" n_reads_unpaired_dropped}" >> input_data.tsv
        echo -e "liftover_fix_pipeline_log\t~{sep="\t" liftover_fix_pipeline_log}" >> input_data.tsv
        echo -e "major_haplogroup\t~{sep="\t" major_haplogroup}" >> input_data.tsv
        echo -e "passes_qc\t~{sep="\t" passes_qc}" >> input_data.tsv
        echo -e "n_liftover_r2_het_ins_sharing_lhs_hom_dele\t~{sep="\t" n_liftover_r2_het_ins_sharing_lhs_hom_dele}" >> input_data.tsv
        echo -e "n_liftover_r2_left_shift\t~{sep="\t" n_liftover_r2_left_shift}" >> input_data.tsv
        echo -e "final_base_level_coverage_metrics\t~{sep="\t" final_base_level_coverage_metrics}" >> input_data.tsv
        echo -e "sex_at_birth\t~{sep="\t" sex_at_birth}" >> input_data.tsv
        echo -e "mt_mean_coverage\t~{sep="\t" mt_mean_coverage}" >> input_data.tsv
        echo -e "r1_vcf_index\t~{sep="\t" r1_vcf_index}" >> input_data.tsv
        echo -e "n_liftover_r2_spanningfixlhs_upstream\t~{sep="\t" n_liftover_r2_spanningfixlhs_upstream}" >> input_data.tsv
        echo -e "nuc_consensus_overlaps\t~{sep="\t" nuc_consensus_overlaps}" >> input_data.tsv
        echo -e "mtdna_consensus_overlaps\t~{sep="\t" mtdna_consensus_overlaps}" >> input_data.tsv
        echo -e "nuc_variants_dropped\t~{sep="\t" nuc_variants_dropped}" >> input_data.tsv
        echo -e "r1_vcf\t~{sep="\t" r1_vcf}" >> input_data.tsv
        echo -e "duplicate_metrics\t~{sep="\t" duplicate_metrics}" >> input_data.tsv
        echo -e "theoretical_sensitivity_metrics\t~{sep="\t" theoretical_sensitivity_metrics}" >> input_data.tsv
        echo -e "cram_path\t~{sep="\t" cram_path}" >> input_data.tsv

        # Process the TSV file in Python
        python3 <<EOF
        import csv

        # Read the input data from the temporary file
        with open("input_data.tsv", "r") as f:
            reader = csv.reader(f, delimiter="\t")
            data = list(reader)

        # Extract headers and rows
        headers = [row[0] for row in data]
        rows = zip(*[row[1:] for row in data])

        # Write the output TSV
        with open("summary.tsv", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(headers)
            writer.writerows(rows)
        EOF
    >>>

    output {
        File summary_tsv = "summary.tsv"
    }

    runtime {
        docker: "python:3.9"
    }
}

task process_tsv_files {
    input {
        File coverage_tsv  # Path to genomics_metrics_Dec142023_1859_02_tz0000.tsv
        File ancestry_tsv  # Path to echo_v4_r2.ancestry_preds.tsv
        File dob_tsv       # Path to echo_DoB_data.tsv
        File input_tsv     # Input TSV file to process
        String output_tsv_name = "processed_data.tsv"  # Name of the output TSV file
        String docker_image = "us.gcr.io/broad-dsde-methods/emeryj/mthail:dockerimageV1"  # Docker image with Python and pandas
    }

    command <<<
        set -euxo pipefail

        python3 <<EOF
        import pandas as pd
        import numpy as np

        # Load the input TSV into a DataFrame
        df = pd.read_csv("~{input_tsv}", sep="\t")

        # Define only the required columns
        columns_needed = ["research_id", "final_base_level_coverage_metrics", "sample", "final_vcf"]

        # Keep only necessary columns
        filtered_df = df

        # Add and rename columns
        filtered_df["s"] = filtered_df["research_id"]
        filtered_df = filtered_df.rename(columns={
        "research_id": "participant_id",
        "final_base_level_coverage_metrics": "coverage",
        "mt_final_vcf": "final_vcf"
        })
        filtered_df_2 = filtered_df.rename(columns={
        "participant_id": "entity:participant_id",
        "mean_coverage": "mt_mean_coverage",
        "median_coverage": "mt_median_coverage",
        "mt_final_vcf": "final_vcf"
        })
        filtered_df_2["freemix_percentage"] = 0

        # Load additional TSV files
        coveragetsv_df = pd.read_csv("~{coverage_tsv}", sep="\t")
        ancestrytsv_df = pd.read_csv("~{ancestry_tsv}", sep="\t")
        dobtsv_df = pd.read_csv("~{dob_tsv}", sep="\t")

        # Merge with filtered_df on 'research_id' and 's'
        filtered_df_2 = filtered_df_2.merge(
        coveragetsv_df[['research_id', 'mean_coverage', 'biosample_collection_date']],
            left_on='s', right_on='research_id', how='left'
        )
        filtered_df_2.drop(columns=['research_id'], inplace=True)
        filtered_df_2 = filtered_df_2.merge(
        ancestrytsv_df[['research_id', 'ancestry_pred']],
            left_on='s', right_on='research_id', how='left'
        )
        filtered_df_2.drop(columns=['research_id'], inplace=True)
        filtered_df_2 = filtered_df_2.merge(
            dobtsv_df[['research_id', 'date_of_birth']],
            left_on='s', right_on='research_id', how='left'
        )
        filtered_df_2.drop(columns=['research_id'], inplace=True)

        # Calculate age
        filtered_df_2['date_of_birth'] = pd.to_datetime(filtered_df_2['date_of_birth'])
        filtered_df_2['biosample_collection_date'] = pd.to_datetime(filtered_df_2['biosample_collection_date'])
        filtered_df_2['age'] = pd.to_numeric(
            np.floor((filtered_df_2['biosample_collection_date'] - filtered_df_2['date_of_birth']).dt.days / 365)
        )
        filtered_df_2['age'] = filtered_df_2['age'].fillna(39).astype(int)

        # Rename columns for compatibility
        filtered_df_2.rename(columns={"mean_coverage": "wgs_mean_coverage"}, inplace=True)
        filtered_df_2.rename(columns={"ancestry_pred": "pop"}, inplace=True)

        # Temporary workaround
        filtered_df_2["wgs_median_coverage"] = filtered_df_2["wgs_mean_coverage"]

        # Filter rows with valid coverage metrics
        filtered_df_2 = filtered_df_2[
            (filtered_df_2['coverage_metrics'].notna()) & (filtered_df_2['coverage_metrics'] != '')
        ]

        # Save the processed DataFrame to a TSV file
        filtered_df_2.to_csv("~{output_tsv_name}", sep="\t", index=False)
        EOF
    >>>

    output {
        File processed_tsv = "~{output_tsv_name}"
    }

    runtime {
        docker: docker_image
        memory: "4 GB"
        cpu: "2"
    }
}

task annotate_coverage {
    input {
        String input_tsv        # Input TSV file path (optional)
        String? output_ht        # Output Hail Table path (optional)
        Int? chunk_size = 100    # Chunk size for combining VCFs (default: 100)
        Boolean overwrite = false  # Overwrite existing files (default: false)
        Boolean keep_targets = false  # Add annotation for target (default: false)
        Int? n_read_partitions = 1  # Number of partitions for reading TSVs (default: 1)
        Boolean hail_only = false  # Skip generating flat files (default: false)
        Int? n_final_partitions = 1000  # Number of partitions for final MT (default: 1000)
        Int? split_merging = 10  # Number of jobs for splitting merging (default: 1)
        String hail_docker      # Docker image with Hail and dependencies
    }

    command <<<
        set -euxo pipefail

        mkdir -p ./tmp
        mkdir -p ./results.ht

        ## TODO - install mtswirl
        ! cd ~;
        ! rm -rf mtSwirl
        ! git clone https://github.com/jamesemery/mtSwirl.git;
        ! cd ~;
        ! pip install --upgrade pip;
        ! pip install --upgrade aiohttp;
        ! pip install --upgrade gnomad;


        # Run the annotate_coverage.py script
        python3 ./mtSwirl/generate_mtdna_call_mt/Terra/annotate_coverage.py \
        ~{if overwrite then "--overwrite" else ""} \
        ~{if keep_targets then "--keep-targets" else ""} \
        --input-tsv ~{input_tsv} \
        --output-ht "./results.ht" \
        --temp-dir "./tmp/" \
        --chunk-size ~{chunk_size} \
        --n-read-partitions ~{n_read_partitions} \
        ~{if hail_only then "--hail-only" else ""} \
        --n-final-partitions ~{n_final_partitions} \
        --split-merging ~{split_merging}

        echo "listing files after running"
        ls -lh ./results.ht*
        tar -czf ./results.ht.tar.gz ./results.ht*
        echo "listing files are tarring"
        ls -lh
        echo "printing directory"
        pwd
        mv ./results.ht.tar.gz /cromwell_root/results.ht.tar.gz
    >>>

    output {
        File output_ht = "results.ht.tar.gz"
    }

    runtime {
        docker: hail_docker
        memory: "8 GB"
        cpu: "24"
        disks: "local-disk 100 SSD"
    }
}