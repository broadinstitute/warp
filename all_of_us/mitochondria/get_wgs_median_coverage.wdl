version 1.0

workflow get_wgs_median_coverage {

    input {
        File wgs_metrics_tsv
        String output_tsv_name = "wgs_median_coverage.tsv"
        Int num_shards = 64
    }

    call split_wgs_tsv as split_input {
        input:
            wgs_metrics_tsv = wgs_metrics_tsv,
            num_shards = num_shards
    }

    scatter (part in split_input.parts) {
        call get_wgs_median_coverage_chunk as chunk_res {
            input:
                wgs_metrics_tsv = part
        }
    }

    call merge_tsvs as merge_results {
        input:
            part_tsvs = chunk_res.median_coverage_tsv,
            output_tsv_name = output_tsv_name
    }

    output {
        File median_coverage_tsv = merge_results.merged_tsv
    }
}

# Split a 2-column TSV (research_id, metrics_path) into N shards, preserving header.
task split_wgs_tsv {
    input {
        File wgs_metrics_tsv
        Int num_shards = 64

        # Runtime parameters
        Int memory_gb = 4
        Int cpu = 1
        Int disk_gb = 20
    }

    command <<<
        set -euxo pipefail

        mkdir -p parts
        HEADER=$(head -n 1 "~{wgs_metrics_tsv}")
        TOTAL_LINES=$(wc -l < "~{wgs_metrics_tsv}")
        if [ "$TOTAL_LINES" -lt 2 ]; then
            echo "Input TSV has no data rows" >&2
            exit 1
        fi

    # Compute data rows and desired shard count
    DATA_LINES=$((TOTAL_LINES - 1))
    NUM=~{num_shards}
    echo "TOTAL_LINES=$TOTAL_LINES DATA_LINES=$DATA_LINES NUM=$NUM" >&2

    # Split data rows into ~NUM balanced chunks using GNU split lines mode
    # This aims for NUM chunks; some may be empty if NUM > DATA_LINES
    tail -n +2 "~{wgs_metrics_tsv}" | split -d -n l/$NUM - parts/shard_
    echo "Produced $(ls parts | wc -l | tr -d ' ') raw parts" >&2
        for f in parts/shard_*; do
            # Prepend header
            mv "$f" "$f.body"
            printf '%s\n' "$HEADER" > "$f.tsv"
            cat "$f.body" >> "$f.tsv"
            rm -f "$f.body"
        done
    >>>

    output {
        Array[File] parts = glob("parts/shard_*.tsv")
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/warp-tools:2.6.1"
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " HDD"
    }
}

# Extract WGS median coverage for a shard of the input TSV.
task get_wgs_median_coverage_chunk {
    input {
        File wgs_metrics_tsv
        String output_tsv_name = "wgs_median_coverage.tsv"

        # Runtime parameters (adjust as needed)
        Int memory_gb = 4
        Int cpu = 1
        Int disk_gb = 20
    }

    command <<<
        set -euxo pipefail

        python3 <<'EOF'
        import sys
        import csv
        import subprocess
        import shlex
        
        in_path = "~{wgs_metrics_tsv}"
        out_path = "~{output_tsv_name}"
        
        # Stream TSV and process two columns per row: research_id, metrics_file_path
        with open(in_path, newline='') as f_in, open(out_path, 'w', newline='') as f_out:
            reader = csv.reader(f_in, delimiter='\t')
            writer = csv.writer(f_out, delimiter='\t')
            header = next(reader, None)
            # Write output header
            writer.writerow(["research_id", "median_coverage"])

            for row in reader:
                if not row:
                    continue
                research_id = row[0]
                metrics_file_path = row[1]
                try:
                    cmd = f"gsutil cat {shlex.quote(metrics_file_path)} | awk '/^[0-9]/ {{ print $4; exit }}'"
                    result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
                    median_coverage = result.stdout.strip()
                except Exception as e:
                    print(f"ERROR: Failed to extract median coverage for {research_id} from {metrics_file_path}: {e}", file=sys.stderr)
                    raise
                writer.writerow([research_id, median_coverage])
        print(f"Wrote {out_path}")
        EOF
        >>>
    
    output {
        File median_coverage_tsv = output_tsv_name
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/warp-tools:2.6.1"
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " HDD"
    }
}

# Merge many per-shard TSVs into a single TSV, keeping a single header.
task merge_tsvs {
    input {
        Array[File] part_tsvs
        String output_tsv_name = "wgs_median_coverage.tsv"

        Int memory_gb = 2
        Int cpu = 1
        Int disk_gb = 10
    }

    command <<<
        set -euxo pipefail

        # Sort file list for reproducibility (avoid executing paths as commands)
        for f in ~{sep=' ' part_tsvs}; do echo "$f"; done | sort > filelist.txt
        first=true
        while IFS= read -r f; do
            if $first; then
                cat "$f" > "~{output_tsv_name}"
                first=false
            else
                tail -n +2 "$f" >> "~{output_tsv_name}"
            fi
        done < filelist.txt
        ls -lh "~{output_tsv_name}"
    >>>

    output {
        File merged_tsv = output_tsv_name
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/warp-tools:2.6.1"
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " HDD"
    }
}
