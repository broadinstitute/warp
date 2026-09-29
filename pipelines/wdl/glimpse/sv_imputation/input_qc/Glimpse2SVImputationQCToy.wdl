version 1.0

struct RuntimeConfig {
    String config_name
    Int cpu
    Int memory_gb
    Int disk_size_gb
    String disk_type
}

workflow ToyInputQC {
    input {
        File gvcf_manifest
        Array[RuntimeConfig] runtime_configs
        String? billing_project_for_rp
    }

    # Extract just the GVCF paths from the manifest to keep the main task clean
    call ParseManifest {
        input:
            gvcf_manifest = gvcf_manifest
    }

    # Run the toy validation task once for each runtime configuration provided
    scatter (config in runtime_configs) {
        call ValidateGvcfsToy {
            input:
                gvcfs = ParseManifest.gvcf_paths,
                runtime_config = config,
                billing_project_for_rp = billing_project_for_rp
        }
    }

    output {
        Array[File] missing_formats_lists = ValidateGvcfsToy.missing_formats
        Array[File] htslib_debug_logs = ValidateGvcfsToy.debug_log
    }
}

task ParseManifest {
    input {
        File gvcf_manifest
    }

    command <<<
        pip install pandas > /dev/null
        python3 -c "import pandas as pd; df = pd.read_csv('~{gvcf_manifest}', sep='\t'); [print(x) for x in df['gvcf_path'].dropna()]" > gvcf_paths.txt
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        cpu: 1
        disks: "local-disk 10 HDD"
        memory: "2 GiB"
        preemptible: 3
    }

    output {
        Array[String] gvcf_paths = read_lines("gvcf_paths.txt")
    }
}

task ValidateGvcfsToy {
    input {
        Array[String] gvcfs
        RuntimeConfig runtime_config
        String? billing_project_for_rp
    }

    String billing_project = select_first([billing_project_for_rp, ""])

    command <<<
        set -uo pipefail
        shopt -s nullglob

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        if [ -n "~{billing_project}" ]; then
            export GCS_REQUESTER_PAYS_PROJECT=~{billing_project}
        fi

        # Maximize htslib and libcurl network verbosity
        export HTS_LOG_LEVEL=trace

        cpu_count=~{runtime_config.cpu}
        printf '%s\n' ~{sep=' ' gvcfs} > all_gvcfs.txt
        mkdir -p chunks results logs
        split -n "r/${cpu_count}" -d --additional-suffix=.txt all_gvcfs.txt chunks/chunk_

        check_gvcf_chunk() {
            local chunk_file="$1"
            local worker_id="$2"
            local gvcfs_with_missing_format_fields=()
            local debug_log="logs/debug_worker_${worker_id}.txt"

            while IFS= read -r gvcf; do
                [ -z "$gvcf" ] && continue
                echo "[worker $worker_id] Validating GVCF file:$gvcf"

                echo -e "\n========================================================" >> "$debug_log"
                echo "STARTING DOWNLOAD: $gvcf" >> "$debug_log"
                echo "TIMESTAMP: $(date -u)" >> "$debug_log"
                
                # Stream the header and capture all trace stderr directly into the debug log
                bcftools view -Ov -h "$gvcf" > "header_${worker_id}.vcf" 2>> "$debug_log"
                local bcf_exit_code=$?
                
                echo "BCFTOOLS EXIT CODE: $bcf_exit_code" >> "$debug_log"

                # Check PL and GT formats
                local format_lines
                format_lines=$(grep '^##FORMAT=<' "header_${worker_id}.vcf" || true)
                local missing_format=false
                
                if ! printf '%s\n' "$format_lines" | grep -q 'ID=PL[,>]'; then
                    missing_format=true
                fi
                if ! printf '%s\n' "$format_lines" | grep -q 'ID=GT[,>]'; then
                    missing_format=true
                fi

                if [ "$missing_format" = true ]; then
                    echo "[worker $worker_id] GVCF file$gvcf is missing PL/GT annotations."
                    gvcfs_with_missing_format_fields+=("$gvcf")
                    echo "RESULT: FAILED PL/GT CHECK" >> "$debug_log"
                else
                    echo "RESULT: PASSED PL/GT CHECK" >> "$debug_log"
                fi
                echo "========================================================" >> "$debug_log"

            done < "$chunk_file"

            if [ ${#gvcfs_with_missing_format_fields[@]} -gt 0 ]; then
                printf '%s\n' "${gvcfs_with_missing_format_fields[@]}" > "results/missing_${worker_id}.txt"
            else
                : > "results/missing_${worker_id}.txt"
            fi
        }

        worker_id=0
        for chunk_file in chunks/chunk_*; do
            check_gvcf_chunk "$chunk_file" "$worker_id" &
            worker_id=$((worker_id + 1))
        done
        wait

        # Aggregate missing formats and debug logs
        cat results/missing_*.txt 2>/dev/null > "~{runtime_config.config_name}_missing_formats.txt" || true
        cat logs/debug_worker_*.txt 2>/dev/null > "~{runtime_config.config_name}_htslib_trace.txt" || true
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/gatk-bcftools-gcloud:1.0.0-4.2.6.1-1.24-1787155398"
        cpu: runtime_config.cpu
        memory: "~{runtime_config.memory_gb} GiB"
        disks: "local-disk ~{runtime_config.disk_size_gb} ~{runtime_config.disk_type}"
        maxRetries: 0
        noAddress: true
    }

    output {
        File missing_formats = "~{runtime_config.config_name}_missing_formats.txt"
        File debug_log = "~{runtime_config.config_name}_htslib_trace.txt"
    }
}
