version 1.0

workflow get_wgs_median_coverage {

    input {
        File wgs_metrics_tsv
        String output_tsv_name = "wgs_median_coverage.tsv"
    }

    call get_wgs_median_coverage as get_median {
        input:
            wgs_metrics_tsv = wgs_metrics_tsv,
            output_tsv_name = output_tsv_name
    }

    output {
        File median_coverage_tsv = get_median.median_coverage_tsv
    }
}

# Extract WGS median coverage from an input TSV.
# NOTE: This task includes a lightweight placeholder implementation.
# You can update the selection/logic once you share the exact column names and rules.
task get_wgs_median_coverage {
    input {
        File wgs_metrics_tsv
        String output_tsv_name = "wgs_median_coverage.tsv"

        # Runtime parameters (adjust as needed)
        Int memory_gb = 8
        Int cpu = 2
        Int disk_gb = 50
    }

    command <<<
        set -euxo pipefail

        python3 <<EOF
        import pandas as pd
        import sys
        import subprocess
        
        in_path = "~{wgs_metrics_tsv}"
        out_path = "~{output_tsv_name}"
        
        # Read the input TSV
        try:
            df = pd.read_csv(in_path, sep='\t')
        except Exception as e:
            print(f"ERROR: Failed to read {in_path}: {e}", file=sys.stderr)
            raise
        
        out_df = pd.DataFrame(columns=['research_id', 'median_coverage'])

        # The df should consist of 2 columns: research_id and a file path specifying the location of the WGS metrics
        # For each row in the df, we want to run the following command using the filepath in the second column: gsutil cat [filepath] | awk '/^[0-9]/ { print $4; exit }'
        
        median_coverages = []
        for index, row in df.iterrows():
            research_id = row[0]
            metrics_file_path = row[1]  # Assuming the second column contains the file path to the WGS metrics         
            try:
                # run the gsutil and awk command to extract the median coverage
                cmd = f"gsutil cat {metrics_file_path} | awk '/^[0-9]/ {{ print $4; exit }}'"
                result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
                median_coverage = result.stdout.strip()
                # write the median_coverage and corresponding research_id to out_df
                out_df.loc[index, 'research_id'] = research_id
                out_df.loc[index, 'median_coverage'] = median_coverage
            except Exception as e:
                print(f"ERROR: Failed to extract median coverage for {research_id} from {metrics_file_path}: {e}", file=sys.stderr)
                raise
        # Write output TSV
        # final output is a
        out_df.to_csv(out_path, sep='\t', index=False)
        print(f"Wrote {out_path} with columns: {list(out_df.columns)}")
        EOF
        >>>
    
    output {
        File median_coverage_tsv = "~{output_tsv_name}"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/python-numpy-pandas:1.0.0-2.2.3-1.25.2"
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " HDD"
    }
}
