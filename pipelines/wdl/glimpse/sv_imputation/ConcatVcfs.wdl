version 1.0

workflow ConcatVcfs {
    # if this changes, update the concat_vcfs_pipeline_version value in Glimpse2SVImputationBatch.wdl
    String pipeline_version = "0.0.1"

    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix
        Boolean do_bcf = true
        Boolean do_sort = false
        String extra_args = "--threads $(nproc) --naive"

        Array[String] regions = []      # if provided, concat within shards and then concat across shards; useful for concat of HiPhase short + SV
        Boolean do_sort_shard = true
        String extra_args_shard = "--threads $(nproc)"
    }

    if (length(regions) > 0) {
        scatter (region in select_first([regions])) {
            call ConcatVcfs as ShardConcatVcfs {
                input:
                    vcfs = vcfs,
                    vcf_idxs = vcf_idxs,
                    output_prefix = output_prefix,
                    do_bcf = do_bcf,
                    do_sort = do_sort_shard,
                    extra_args = extra_args_shard,
                    region = region
            }
        }
    }

    call ConcatVcfs {
        input:
            vcfs = select_first([ShardConcatVcfs.concatenated_vcf, vcfs]),
            vcf_idxs = select_first([ShardConcatVcfs.concatenated_vcf_idx, vcf_idxs]),
            output_prefix = output_prefix,
            do_bcf = do_bcf,
            do_sort = do_sort,
            extra_args = extra_args
    }

    output {
        File concatenated_vcf = ConcatVcfs.concatenated_vcf
        File concatenated_vcf_idx = ConcatVcfs.concatenated_vcf_idx
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

task ConcatVcfs {
    input{
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix
        Boolean do_bcf = true
        Boolean do_sort = false
        String? extra_args
        String? region

        RuntimeAttr? runtime_attr_override
    }

    # If sorting, provide extra disk space for the temporary sort shards
    Int disk_gb = if do_sort then 10 + 4 * ceil(size(vcfs, "GiB")) else 10 + 2 * ceil(size(vcfs, "GiB"))

    String format_arg = if do_bcf then "-Ob" else "-Oz"
    String suffix = if do_bcf then "bcf" else "vcf.gz"
    String idx_suffix = if do_bcf then "bcf.csi" else "vcf.gz.tbi"
    String index_arg = if do_bcf then "-c" else "-t"

    command <<<
        set -euox pipefail

        # Start zero-overhead background heartbeat monitor
        (
            echo "Starting concat monitoring..." >&2
            while true; do
                if [ -f "~{output_prefix}.~{suffix}" ]; then
                    # Phase 2/No-Sort: Final file is being written
                    SIZE=$(ls -lh "~{output_prefix}.~{suffix}" | awk '{print $5}')
                    echo "[Heartbeat] Final output ~{output_prefix}.~{suffix} is currently $SIZE..." >&2
                elif [ -d "sort_tmp_dir" ]; then
                    # Phase 1 (Sorting): Temp directory is filling up
                    SIZE=$(du -sh sort_tmp_dir | awk '{print $1}')
                    echo "[Heartbeat] Sorting in progress. Temp shards total $SIZE..." >&2
                else
                    echo "[Heartbeat] Processing started, waiting for I/O..." >&2
                fi
                sleep 60
            done
        ) &
        HEARTBEAT_PID=$!

        if [ "~{do_sort}" == "true" ]; then
            echo "Concatenating and piping to bcftools sort..."

            # Create a dedicated temp directory so we can monitor its size
            mkdir sort_tmp_dir

            # Output as uncompressed BCF (-Ou) and route sort temp files to our directory (-T)
            bcftools concat ~{"--regions-overlap 0 -r " + region} \
                -f ~{write_lines(vcfs)} \
                ~{extra_args} \
                -Ou | bcftools sort -m 2G -T sort_tmp_dir/tmp ~{format_arg} -o ~{output_prefix}.~{suffix}
        else
            echo "Concatenating directly to disk..."
            bcftools concat ~{"--regions-overlap 0 -r " + region} \
                -f ~{write_lines(vcfs)} \
                ~{extra_args} \
                ~{format_arg} -o ~{output_prefix}.~{suffix}
        fi

        bcftools index ~{index_arg} ~{output_prefix}.~{suffix}

        # Kill the background monitor the second the pipeline finishes
        kill $HEARTBEAT_PID || true
    >>>

    output {
        File concatenated_vcf = "~{output_prefix}.~{suffix}"
        File concatenated_vcf_idx = "~{output_prefix}.~{idx_suffix}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
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
