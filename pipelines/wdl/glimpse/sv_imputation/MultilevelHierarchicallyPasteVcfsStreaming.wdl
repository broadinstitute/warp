version 1.0

# NOTE: We assume we are merging squared-off single-sample VCFs (enforced by checking that the number of records in each VCF is identical when localizing shards with bcftools view)

workflow MultilevelHierarchicallyMergeVcfs {
    # if this changes, update the multi_level_paste_pipeline_version value in PreprocessPLsGVCF.wdl
    String pipeline_version = "0.0.1"

    input {
        Array[String]? vcfs_array
        Array[String]? vcf_idxs_array
        File? vcfs_fofn
        File? vcf_idxs_fofn
        Array[String] regions   # bcftools regions, e.g. ["chr1,chr2,chr3", "chr4,chr5,chr6", ...]
        Array[Int] batch_sizes  # Parameterizable hierarchical levels, e.g., [100, 50]
        Array[Boolean] do_localization # Whether to localize at each corresponding level
        Array[Int] timeouts_min  # Timeouts in minutes per level. Set to 0 to disable. e.g., [720, 720]
        String output_prefix

        String extra_merge_args = "--threads $(nproc) --info ID,RAF --format GT,DS,GP"

        String extra_concat_args = "--threads $(nproc) --naive"
    }

    Array[String] vcfs_in = if defined(vcfs_array) then select_first([vcfs_array]) else read_lines(select_first([vcfs_fofn]))
    Array[String] vcf_idxs_in = if defined(vcf_idxs_array) then select_first([vcf_idxs_array]) else read_lines(select_first([vcf_idxs_fofn]))

    call CreateBatches as L0_Batches {
        input:
            vcfs = vcfs_in,
            vcf_idxs = vcf_idxs_in,
            batch_size = batch_sizes[0]
    }

    # Scatter by region FIRST to isolate chunks and reduce combinatorial explosion
    scatter (j in range(length(regions))) {
        String region = regions[j]
        String region_prefix = output_prefix + ".region-" + j

        # ==========================================
        # LEVEL 0
        # ==========================================
        scatter (i in range(length(L0_Batches.vcf_batch_fofns))) {
            call MergeVcfs as L0_Merge {
                input:
                    vcfs_localize = if do_localization[0] then read_lines(L0_Batches.vcf_batch_fofns[i]) else [],
                    vcf_idxs_localize = if do_localization[0] then read_lines(L0_Batches.vcf_idx_batch_fofns[i]) else [],
                    vcfs_stream = if !do_localization[0] then read_lines(L0_Batches.vcf_batch_fofns[i]) else [],
                    vcf_idxs_stream = if !do_localization[0] then read_lines(L0_Batches.vcf_idx_batch_fofns[i]) else [],
                    timeout_min = timeouts_min[0],
                    region = region,
                    output_prefix = region_prefix + ".L0-" + i,
                    extra_args = "-r " + region + " " + extra_merge_args
            }
        }

        Array[File] l0_vcfs = L0_Merge.merged_vcf
        Array[File] l0_idxs = L0_Merge.merged_vcf_idx

        # ==========================================
        # LEVEL 1
        # ==========================================
        if (length(batch_sizes) > 1) {
            call CreateBatches as L1_Batches {
                input:
                    vcfs = l0_vcfs,
                    vcf_idxs = l0_idxs,
                    batch_size = batch_sizes[1]
            }

            scatter (i in range(length(L1_Batches.vcf_batch_fofns))) {
                call MergeVcfs as L1_Merge {
                    input:
                        vcfs_localize = if do_localization[1] then read_lines(L1_Batches.vcf_batch_fofns[i]) else [],
                        vcf_idxs_localize = if do_localization[1] then read_lines(L1_Batches.vcf_idx_batch_fofns[i]) else [],
                        vcfs_stream = if !do_localization[1] then read_lines(L1_Batches.vcf_batch_fofns[i]) else [],
                        vcf_idxs_stream = if !do_localization[1] then read_lines(L1_Batches.vcf_idx_batch_fofns[i]) else [],
                        timeout_min = timeouts_min[1],
                        region = region,
                        output_prefix = region_prefix + ".L1-" + i,
                        extra_args = "-r " + region + " " + extra_merge_args
                }
            }
        }

        Array[File] l1_vcfs = select_first([L1_Merge.merged_vcf, l0_vcfs])
        Array[File] l1_idxs = select_first([L1_Merge.merged_vcf_idx, l0_idxs])

        # ==========================================
        # LEVEL 2
        # ==========================================
        if (length(batch_sizes) > 2) {
            call CreateBatches as L2_Batches {
                input:
                    vcfs = l1_vcfs,
                    vcf_idxs = l1_idxs,
                    batch_size = batch_sizes[2]
            }

            scatter (i in range(length(L2_Batches.vcf_batch_fofns))) {
                call MergeVcfs as L2_Merge {
                    input:
                        vcfs_localize = if do_localization[2] then read_lines(L2_Batches.vcf_batch_fofns[i]) else [],
                        vcf_idxs_localize = if do_localization[2] then read_lines(L2_Batches.vcf_idx_batch_fofns[i]) else [],
                        vcfs_stream = if !do_localization[2] then read_lines(L2_Batches.vcf_batch_fofns[i]) else [],
                        vcf_idxs_stream = if !do_localization[2] then read_lines(L2_Batches.vcf_idx_batch_fofns[i]) else [],
                        timeout_min = timeouts_min[2],
                        region = region,
                        output_prefix = region_prefix + ".L2-" + i,
                        extra_args = "-r " + region + " " + extra_merge_args
                }
            }
        }

        Array[File] l2_vcfs = select_first([L2_Merge.merged_vcf, l1_vcfs])
        Array[File] l2_idxs = select_first([L2_Merge.merged_vcf_idx, l1_idxs])

        # ==========================================
        # FINAL REGION COLLAPSE
        # ==========================================
        # Safety Net: Defaults to localizing workspace intermediate files, no timeout applied (0).
        if (length(l2_vcfs) > 1) {
            call MergeVcfs as FinalRegionMerge {
                input:
                    vcfs_localize = l2_vcfs,
                    vcf_idxs_localize = l2_idxs,
                    timeout_min = 0,
                    region = region,
                    output_prefix = region_prefix + ".final",
                    extra_args = "-r " + region + " " + extra_merge_args
            }
        }

        # Select exactly 1 file for this region to pass to the final concat step
        File final_region_vcf = select_first([FinalRegionMerge.merged_vcf, l2_vcfs[0]])
        File final_region_idx = select_first([FinalRegionMerge.merged_vcf_idx, l2_idxs[0]])
    }

    # concatenate all regions together
    call ConcatVcfs {
        input:
            vcfs = final_region_vcf,
            vcf_idxs = final_region_idx,
            output_prefix = output_prefix,
            extra_args = extra_concat_args
    }

    output {
        File merged_vcf = ConcatVcfs.concatenated_vcf
        File merged_vcf_idx = ConcatVcfs.concatenated_vcf_idx
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    String? disk_type
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

task CreateBatches {
    input {
        Array[String] vcfs
        Array[String] vcf_idxs
        Int batch_size

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euox pipefail

        # Split with -a 4 guarantees strict alphanumeric ordering up to 456,976 batches
        cat ~{write_lines(vcfs)} | split -a 4 -l ~{batch_size} - vcf_batch_
        cat ~{write_lines(vcf_idxs)} | split -a 4 -l ~{batch_size} - vcf_idx_batch_
    >>>

    output {
        Array[File] vcf_batch_fofns = glob("vcf_batch_*")
        Array[File] vcf_idx_batch_fofns = glob("vcf_idx_batch_*")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            10,
        boot_disk_gb:       10,
        disk_type:          "HDD",
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/bcftools-vcftools:sps_sv_docker_images"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MergeVcfs {
    input {
        Array[File] vcfs_localize = []
        Array[File] vcf_idxs_localize = []
        Array[String] vcfs_stream = []
        Array[String] vcf_idxs_stream = []

        Int timeout_min
        String? region
        String output_prefix
        String? extra_args

        RuntimeAttr? runtime_attr_override
    }

    # Dynamically sizes disk if localizing, defaults to 50GB if streaming
    Int disk_gb = if length(vcfs_localize) > 0 then 10 + 2 * ceil(size(vcfs_localize, "GiB")) else 50

    command <<<
        set -euox pipefail

        if [ ~{length(vcfs_localize)} -gt 0 ]; then
            echo "Localizing files natively via Cromwell..."
            cat ~{write_lines(vcfs_localize)} > merge_list.txt
        else
            echo "Slice-and-downloading regions via bcftools view..."
            export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

            mkdir -p localized_subsets

            # Stitch the VCF and IDX strings together using htslib's explicit index syntax
            paste \
                ~{write_lines(vcfs_stream)} \
                ~{write_lines(vcf_idxs_stream)} \
                | awk '{print $1"##idx##"$2}' > remote_list.txt

            # Prepend the line number (NR) using awk, separated by a pipe, and pass to xargs
            awk '{print NR"|"$0}' remote_list.txt | xargs -P $(nproc) -I {} bash -c '
                set -euox pipefail

                # Split the line number and the URL
                IFS="|" read -r LINE_NUM URL <<< "{}"

                # Extract the base filename, ignoring the index syntax
                BASENAME=$(basename "${URL%%##idx##*}")
                OUT_BCF="localized_subsets/${LINE_NUM}_${BASENAME}"

                # Setup timeout prefix if timeout_min > 0
                TIMEOUT_CMD=""
                if [ ~{timeout_min} -gt 0 ]; then
                    TIMEOUT_CMD="timeout ~{timeout_min}m"
                fi

                # Retry loop to catch transient network hangs
                MAX_RETRIES=3
                for i in $(seq 1 $MAX_RETRIES); do
                    echo "Attempt $i for ${BASENAME}..."

                    # Wrap in timeout and use --write-index to generate the .csi simultaneously
                    if ${TIMEOUT_CMD} bcftools view \
                        ~{if defined(region) then "--regions-overlap 0 -r " + region else ""} \
                        --write-index \
                        "${URL}" -Ob -o "${OUT_BCF}"; then

                        exit 0
                    else
                        echo "WARNING: Download failed or timed out for ${URL}. Retrying in 5 seconds..." >&2
                        sleep 5
                    fi
                done

                echo "FATAL: Failed to download ${URL} after $MAX_RETRIES attempts." >&2
                exit 1
            '

            # Build the merge list matching the exact 1-based line number (NR) prefix we just created
            awk -F '##idx##' '{
                n = split($1, a, "/");
                print "localized_subsets/" NR "_" a[n]
            }' remote_list.txt > merge_list.txt

            # ==========================================
            # FAST RECORD COUNT ASSERTION
            # ==========================================
            echo "Verifying record counts via .csi indices..."
            EXPECTED_RECORDS=""

            while IFS= read -r file; do
                # Dynamically extract record counts from either a CSI or a Tabix index mapping safely
                if [[ "$file" == *.vcf.gz ]]; then
                    RECORDS=$(bcftools index -n "$file" 2>/dev/null || tabix -l "$file" | wc -l)
                else
                    RECORDS=$(bcftools index -n "$file")
                fi

                if [ -z "$EXPECTED_RECORDS" ]; then
                    EXPECTED_RECORDS=$RECORDS
                    echo "Baseline established: Expecting exactly $EXPECTED_RECORDS records per file."
                elif [ "$RECORDS" != "$EXPECTED_RECORDS" ]; then
                    echo "FATAL: Record count mismatch! $file has $RECORDS records, expected $EXPECTED_RECORDS." >&2
                    exit 1
                fi
            done < merge_list.txt

            echo "SUCCESS: All localized subsets perfectly match at $EXPECTED_RECORDS records."
        fi

        # Start a zero-overhead background heartbeat monitor
        (
            echo "Starting merge monitoring..." >&2
            while true; do
                if [ -f "~{output_prefix}.bcf" ]; then
                    # Fetch the human-readable file size safely
                    SIZE=$(ls -lh "~{output_prefix}.bcf" | awk '{print $5}')
                    echo "[Heartbeat] ~{output_prefix}.bcf is currently $SIZE..." >&2
                fi
                sleep 60
            done
        ) &
        HEARTBEAT_PID=$!

        # ==========================================
        # EXECUTE CUSTOM MERGE
        # ==========================================
        # Execute the compiled tool, pasting positional inputs straight from our list
        /usr/local/bin/paste-vcfs \
            ~{extra_args} \
            -o ~{output_prefix}.bcf \
            $(cat merge_list.txt)

        bcftools index ~{output_prefix}.bcf

        kill $HEARTBEAT_PID || true
    >>>

    output {
        File merged_vcf = "~{output_prefix}.bcf"
        File merged_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsde-methods/sshah/sv_rust_tools:1.0.0-5dc0f19-1784242277"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ConcatVcfs {
    input{
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix
        String? extra_args

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size(vcfs, "GiB"))

    command <<<
        set -euox pipefail

        bcftools concat \
            -f ~{write_lines(vcfs)} \
            ~{extra_args} \
            -Ob -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File concatenated_vcf = "~{output_prefix}.bcf"
        File concatenated_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-gotc-prod/bcftools-vcftools:sps_sv_docker_images"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
