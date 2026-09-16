version 1.0

task validate_rnaseqc_manifests {
    input {
        File sample_manifest
        String prefix
        Int batch_size

        String docker_image
        Int memory_gb
        Int disk_space_gb
        Int num_preempt
    }

    # Keep file creation in task scope for Terra's backend filesystem.
    # Serialize the string without inserting its contents into shell source.
    File prefix_input_file = write_lines([prefix])

    command <<<
        set -euo pipefail

        log() {
            printf '[%s] %s\n' "$(date '+%b %d %H:%M:%S')" "$1"
        }

        log "stage=validate_manifest start_time=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        log "Validating the combined RNA-SeQC sample manifest"
        python3 /opt/prepare_qtl/scripts/expression/rnaseqc/merge_rnaseqc.py validate-manifest \
            --input "~{sample_manifest}" \
            --batch-size ~{batch_size} \
            --prefix-file "~{prefix_input_file}"

        log "Saving the validated output prefix"
        cp -- "~{prefix_input_file}" prefix.txt

        sample_count=$(<sample_count.txt)
        batch_count=$(<batch_count.txt)
        log "dimensions=samples:$sample_count,batches:$batch_count"
        log "outputs=sample_count.txt,batch_count.txt,include_insert_sizes.txt,prefix.txt"
        log "stage=validate_manifest completion_time=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    >>>

    output {
        Int sample_count = read_int("sample_count.txt")
        Int batch_count = read_int("batch_count.txt")
        Boolean include_insert_sizes = read_boolean("include_insert_sizes.txt")
        File prefix_file = "prefix.txt"
    }

    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_space_gb} HDD"
        cpu: 1
        preemptible: num_preempt
    }
}

task aggregate_rnaseqc_batch {
    input {
        File sample_manifest
        File prefix_file
        Int batch_index
        Int batch_size
        Boolean include_insert_sizes
        Boolean merge_exons

        String docker_image
        Int memory_gb
        Int disk_space_gb
        Int num_threads
        Int num_preempt
    }

    command <<<
        set -euo pipefail
        export LC_ALL=C

        log() {
            printf '[%s] %s\n' "$(date '+%b %d %H:%M:%S')" "$1"
        }

        batch_number=$((~{batch_index} + 1))
        root_prefix=$(<"~{prefix_file}")
        batch_prefix="${root_prefix}.batch_~{batch_index}"

        if [[ ~{num_threads} -lt 1 ]]; then
            printf 'ERROR: num_threads must be at least 1\n' >&2
            exit 1
        fi

        log "stage=aggregate_batch batch=$batch_number start_time=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        log "merge_exons=~{merge_exons}"
        log "Preparing batch $batch_number"
        mkdir -p individual_outputs
        python3 /opt/prepare_qtl/scripts/expression/rnaseqc/merge_rnaseqc.py prepare-batch \
            --input "~{sample_manifest}" \
            --batch-index ~{batch_index} \
            --batch-size ~{batch_size} \
            --staging-directory individual_outputs~{if merge_exons then "" else " --skip-exons"}

        batch_sample_count=$(awk 'END { print NR + 0 }' batch_sample_ids.list)
        log "dimensions=samples:$batch_sample_count,batch:$batch_number"

        log "Staging $batch_sample_count samples for batch $batch_number"
        # shellcheck disable=SC2016
        while IFS=$'\t' read -r source destination; do
            printf '%s\0%s\0' "$source" "$destination"
        done < transfers.tsv | xargs -0 -n 2 -P ~{num_threads} bash -c 'gsutil cp "$1" "$2"' _

        log "Merging TPM GCT files for batch $batch_number"
        python3 /opt/prepare_qtl/scripts/expression/rnaseqc/merge_rnaseqc.py gct \
            --input-list local_tpm.list \
            --output "$batch_prefix.gene_tpm.gct.gz" \
            --sample-output batch_samples.txt \
            --sample-names batch_sample_ids.list

        log "Merging count GCT files for batch $batch_number"
        python3 /opt/prepare_qtl/scripts/expression/rnaseqc/merge_rnaseqc.py gct \
            --input-list local_count.list \
            --output "$batch_prefix.gene_reads.gct.gz" \
            --sample-names batch_sample_ids.list

        if [[ "~{merge_exons}" == "true" ]]; then
            log "Merging exon-count GCT files for batch $batch_number"
            python3 /opt/prepare_qtl/scripts/expression/rnaseqc/merge_rnaseqc.py gct \
                --input-list local_exon.list \
                --output "$batch_prefix.exon_reads.gct.gz" \
                --sample-names batch_sample_ids.list
        else
            log "Skipping exon-count downloads and merging for batch $batch_number"
        fi

        log "Merging metrics files for batch $batch_number"
        python3 /opt/prepare_qtl/scripts/expression/rnaseqc/merge_rnaseqc.py metrics-individual \
            --input-list local_metrics.list \
            --output "$batch_prefix.metrics.txt.gz" \
            --sample-names batch_sample_ids.list

        if [[ "~{include_insert_sizes}" == "true" ]]; then
            log "Merging insert-size files for batch $batch_number"
            python3 /opt/prepare_qtl/scripts/expression/rnaseqc/merge_rnaseqc.py insert-sizes-individual \
                --input-list local_insert.list \
                --expected-samples batch_samples.txt \
                --output "$batch_prefix.insert_size_hists.txt.gz"
        fi

        log "outputs=$batch_prefix.gene_tpm.gct.gz,$batch_prefix.gene_reads.gct.gz,$batch_prefix.metrics.txt.gz"
        if [[ "~{merge_exons}" == "true" ]]; then
            log "outputs=$batch_prefix.exon_reads.gct.gz"
        fi
        if [[ "~{include_insert_sizes}" == "true" ]]; then
            log "outputs=$batch_prefix.insert_size_hists.txt.gz"
        fi
        log "stage=aggregate_batch batch=$batch_number completion_time=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    >>>

    output {
        File metrics = glob("*.metrics.txt.gz")[0]
        File tpm_gct = glob("*.gene_tpm.gct.gz")[0]
        File count_gct = glob("*.gene_reads.gct.gz")[0]
        Array[File] exon_count_gct = glob("*.exon_reads.gct.gz")
        Array[File] insert_size_hists = glob("*.insert_size_hists.txt.gz")
    }

    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_space_gb} HDD"
        cpu: num_threads
        preemptible: num_preempt
    }
}

task merge_rnaseqc_batches {
    input {
        Array[File] batch_tpm_gcts
        Array[File] batch_count_gcts
        Array[File] batch_exon_count_gcts
        Array[File] batch_metrics
        Array[File] batch_insert_size_hists
        File prefix_file
        Boolean include_insert_sizes
        Boolean merge_exons

        String docker_image
        Int memory_gb
        Int disk_space_gb
        Int num_preempt
    }

    command <<<
        set -euo pipefail
        export LC_ALL=C

        log() {
            printf '[%s] %s\n' "$(date '+%b %d %H:%M:%S')" "$1"
        }

        prefix=$(<"~{prefix_file}")

        log "stage=merge_cohort start_time=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        log "merge_exons=~{merge_exons}"
        log "Merging batch-level TPM GCT files"
        python3 /opt/prepare_qtl/scripts/expression/rnaseqc/merge_rnaseqc.py gct \
            --input-list "~{write_lines(batch_tpm_gcts)}" \
            --output "$prefix.gene_tpm.gct.gz" \
            --sample-output cohort_samples.txt

        sample_count=$(awk 'END { print NR + 0 }' cohort_samples.txt)
        log "dimensions=samples:$sample_count,batches:~{length(batch_tpm_gcts)}"

        log "Merging batch-level count GCT files"
        python3 /opt/prepare_qtl/scripts/expression/rnaseqc/merge_rnaseqc.py gct \
            --input-list "~{write_lines(batch_count_gcts)}" \
            --expected-samples cohort_samples.txt \
            --output "$prefix.gene_reads.gct.gz"

        if [[ "~{merge_exons}" == "true" ]]; then
            log "Merging batch-level exon-count GCT files"
            python3 /opt/prepare_qtl/scripts/expression/rnaseqc/merge_rnaseqc.py gct \
                --input-list "~{write_lines(batch_exon_count_gcts)}" \
                --expected-samples cohort_samples.txt \
                --output "$prefix.exon_reads.gct.gz"
        else
            log "Skipping the final exon-count merge"
        fi

        log "Merging batch-level metrics files"
        python3 /opt/prepare_qtl/scripts/expression/rnaseqc/merge_rnaseqc.py metrics-aggregated \
            --input-list "~{write_lines(batch_metrics)}" \
            --expected-samples cohort_samples.txt \
            --output "$prefix.metrics.txt.gz"

        if [[ "~{include_insert_sizes}" == "true" ]]; then
            log "Merging batch-level insert-size files"
            python3 /opt/prepare_qtl/scripts/expression/rnaseqc/merge_rnaseqc.py insert-sizes-aggregated \
                --input-list "~{write_lines(batch_insert_size_hists)}" \
                --expected-samples cohort_samples.txt \
                --output "$prefix.insert_size_hists.txt.gz"
        fi

        log "outputs=$prefix.gene_tpm.gct.gz,$prefix.gene_reads.gct.gz,$prefix.metrics.txt.gz"
        if [[ "~{merge_exons}" == "true" ]]; then
            log "outputs=$prefix.exon_reads.gct.gz"
        fi
        if [[ "~{include_insert_sizes}" == "true" ]]; then
            log "outputs=$prefix.insert_size_hists.txt.gz"
        fi
        log "stage=merge_cohort completion_time=$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    >>>

    output {
        File metrics = glob("*.metrics.txt.gz")[0]
        File tpm_gct = glob("*.gene_tpm.gct.gz")[0]
        File count_gct = glob("*.gene_reads.gct.gz")[0]
        Array[File] exon_count_gct = glob("*.exon_reads.gct.gz")
        Array[File] insert_size_hists = glob("*.insert_size_hists.txt.gz")
    }

    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_space_gb} HDD"
        cpu: 1
        preemptible: num_preempt
    }
}

workflow rnaseqc2_aggregate_batched_workflow {
    input {
        File sample_manifest
        String prefix
        Boolean merge_exons

        Int batch_size = 100
        String docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-rnaseqc2-aggregation@sha256:cb725753b77ff558d5390966c0faea4ada431845f36239b00e9b3e2a012faa1f"
        Int validation_memory_gb = 1
        Int validation_disk_space_gb = 10
        Int batch_memory_gb = 4
        Int batch_disk_space_gb = 100
        Int merge_memory_gb = 8
        Int merge_disk_space_gb
        Int num_threads = 8
        Int num_preempt = 2
    }
    String pipeline_version = "aou_9.1.0"

    call validate_rnaseqc_manifests {
        input:
            sample_manifest = sample_manifest,
            prefix = prefix,
            batch_size = batch_size,
            docker_image = docker_image,
            memory_gb = validation_memory_gb,
            disk_space_gb = validation_disk_space_gb,
            num_preempt = num_preempt
    }

    scatter (batch_index in range(validate_rnaseqc_manifests.batch_count)) {
        call aggregate_rnaseqc_batch {
            input:
                sample_manifest = sample_manifest,
                prefix_file = validate_rnaseqc_manifests.prefix_file,
                batch_index = batch_index,
                batch_size = batch_size,
                include_insert_sizes = validate_rnaseqc_manifests.include_insert_sizes,
                merge_exons = merge_exons,
                docker_image = docker_image,
                memory_gb = batch_memory_gb,
                disk_space_gb = batch_disk_space_gb,
                num_threads = num_threads,
                num_preempt = num_preempt
        }
    }

    call merge_rnaseqc_batches {
        input:
            batch_tpm_gcts = aggregate_rnaseqc_batch.tpm_gct,
            batch_count_gcts = aggregate_rnaseqc_batch.count_gct,
            batch_exon_count_gcts = flatten(aggregate_rnaseqc_batch.exon_count_gct),
            batch_metrics = aggregate_rnaseqc_batch.metrics,
            batch_insert_size_hists = flatten(aggregate_rnaseqc_batch.insert_size_hists),
            prefix_file = validate_rnaseqc_manifests.prefix_file,
            include_insert_sizes = validate_rnaseqc_manifests.include_insert_sizes,
            merge_exons = merge_exons,
            docker_image = docker_image,
            memory_gb = merge_memory_gb,
            disk_space_gb = merge_disk_space_gb,
            num_preempt = num_preempt
    }

    output {
        File metrics = merge_rnaseqc_batches.metrics
        File tpm_gct = merge_rnaseqc_batches.tpm_gct
        File count_gct = merge_rnaseqc_batches.count_gct
        Array[File] exon_count_gct = merge_rnaseqc_batches.exon_count_gct
        Array[File] insert_size_hists = merge_rnaseqc_batches.insert_size_hists
        Int sample_count = validate_rnaseqc_manifests.sample_count
        Int batch_count = validate_rnaseqc_manifests.batch_count
    }

    meta {
        description: "Aggregate RNA-SeQC outputs in bounded batches and stream the final cohort merge."
    }
}
