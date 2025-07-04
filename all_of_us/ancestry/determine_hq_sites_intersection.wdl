version 1.0
# Determine the high-quality sites as an intersection of the training high quality sites (from determine hq_sites.wdl) and
#  the input data (ordered_vcf_shards).
workflow determine_hq_sites_intersection {
    input {

        # The high quality sites, with the samples.  These are generated by select_samples_for_training.wdl.
        #  I.e. the high-quality sites w/ the training samples GT
        File training_vcf_bgz
        File training_vcf_bgz_idx

        # The high quality sites, without the samples (i.e. sites-only).  These are generated by select_samples_for_training.wdl.
        #  This helps with speed of this workflow.  Must correspond to the training_vcf_bgz above
        File training_vcf_so_bgz
        File training_vcf_so_bgz_idx

        # The data we wish to obtain ancestry
        # These should have the same samples in each shard, but consecutive regions of the genome
        Array[File] ordered_vcf_shards_in
        Array[File] ordered_vcf_shards_idx_in

        # Optional file that specified file of file names.  If specified, the ordered_vcf_shards_in parameter is ignored
        File? ordered_vcf_shards_list
        File? ordered_vcf_shards_idx_list

        String final_output_prefix

        # MUST be a valid GS URL
        String? service_account_json

        # Further intersect sites.  Eg, exome intervals.  This can be used to further reduce the area to search for SNPs
        File? intersecting_intervals
    }

    Array[File] ordered_vcf_shards = if (defined(ordered_vcf_shards_list)) then read_lines(select_first([ordered_vcf_shards_list, ""])) else ordered_vcf_shards_in
    Array[File] ordered_vcf_shards_idx = if (defined(ordered_vcf_shards_idx_list)) then read_lines(select_first([ordered_vcf_shards_idx_list, ""])) else ordered_vcf_shards_idx_in
    String pipeline_version = "aou-8.0.0"

    # Get the high quality sites that are called in the test data (intersection file).
    #  Return as a full VCF of the training data ("hq full").  And return the count as well.
    Array[Pair[File,File]] vcf_shard_pairs = zip(ordered_vcf_shards, ordered_vcf_shards_idx)
    scatter (i in range(length(vcf_shard_pairs))) {
        Pair[File,File] vcf_shard_pair = vcf_shard_pairs[i]

        # Note that we do not do LD pruning on the VCF shards.
        call sitesOnlyAndHQFilterVcf as createDataSites {
            input:
                vcf = vcf_shard_pair.left,
                vcf_idx = vcf_shard_pair.right,
                service_account_json=service_account_json,
                vcf_intervals = training_vcf_so_bgz,
                vcf_intervals_idx = training_vcf_so_bgz_idx,
                intersecting_intervals=intersecting_intervals,
                id = i
        }

        call intersect_vcfs_as_sites_only {
            input:
                vcf1 = createDataSites.vcf_sites_only,
                vcf1_index = createDataSites.vcf_sites_only_idx,
                vcf2 = training_vcf_so_bgz,
                vcf2_index = training_vcf_so_bgz_idx,
                output1_basename = "sites_only_intersection." + i
        }
    }

    call merge_vcf_bgzs as merge_sites_only_intersection {
        input:
            input_vcfs = intersect_vcfs_as_sites_only.intersected_vcf,
            input_vcf_indices = intersect_vcfs_as_sites_only.intersected_vcf_idx,
            output_name = "merged_sites_only_intersection"
    }

    call filter_by_sites_only as filter_training_set {
        input:
            vcf=training_vcf_bgz,
            vcf_idx=training_vcf_bgz_idx,
            sites_only_vcf = merge_sites_only_intersection.merged_vcf,
            sites_only_vcf_idx = merge_sites_only_intersection.merged_vcf_idx,
            output_name="full_training_sites_filtered",
            id="0"
    }


    scatter (j in range(length(vcf_shard_pairs))) {
        Pair[File,File] vcf_shard_pair2 = vcf_shard_pairs[j]
        call filter_by_sites_only as filter_vcf_shard {
            input:
                vcf=vcf_shard_pair2.left,
                vcf_idx=vcf_shard_pair2.right,
                sites_only_vcf = merge_sites_only_intersection.merged_vcf,
                sites_only_vcf_idx = merge_sites_only_intersection.merged_vcf_idx,
                output_name="full_data_sites_filtered",
                service_account_json=service_account_json,
                id = j
        }
    }

    call merge_vcf_bgzs as merge_vcf_shards {
        input:
            input_vcfs = filter_vcf_shard.filtered_vcf,
            input_vcf_indices = filter_vcf_shard.filtered_vcf_idx,
            output_name = "merged_data_shards"
    }

    output {
        # A sites-only file of the variants that are the training hq sites minus sites that did not appear in the data
        File hq_variants_intersection = merge_sites_only_intersection.merged_vcf
        File hq_variants_intersection_idx = merge_sites_only_intersection.merged_vcf_idx

        # The merged VCF of data files at the hq sites (i.e. at the hq_variants_intersection)
        File merged_vcf_shards = merge_vcf_shards.merged_vcf
        File merged_vcf_shards_idx = merge_vcf_shards.merged_vcf_idx

        # The VCF of training samples at the hq sites (i.e. at the hq_variants_intersection)
        File filtered_training_set = filter_training_set.filtered_vcf
        File filtered_training_set_idx = filter_training_set.filtered_vcf_idx
    }
}

task sitesOnlyAndHQFilterVcf {
    input {
        File vcf
        File vcf_idx
        File vcf_intervals
        File vcf_intervals_idx

        # Any String can go here, but we recommend the shard index if possible.
        String id
        String? service_account_json

        # Further intersect sites.  Eg, exome intervals
        File? intersecting_intervals
    }
    String output_filename = basename(vcf) + "." + id + ".sites_only.vcf.gz"
    String output_filename_idx = output_filename + ".tbi"
    String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'
    String service_account_basename_pre = if (defined(service_account_json)) then service_account_json else ''
    String service_account_basename = basename(service_account_basename_pre)
    String input_vcf_basename = basename(vcf)
    String updated_input_vcf = if (defined(service_account_json)) then input_vcf_basename else vcf

    parameter_meta {
        vcf: {localization_optional: true}
        vcf_idx: {localization_optional: true}
    }

    command <<<
        set -e

        if [ ~{has_service_account_file} = 'true' ]; then
        gsutil cp ~{service_account_json} ~{service_account_basename}
        export GOOGLE_APPLICATION_CREDENTIALS=~{service_account_basename}
        gcloud auth activate-service-account --key-file='~{service_account_basename}'
        gsutil cp ~{vcf} .
        gsutil cp ~{vcf_idx} .
        fi

        gatk --java-options "-Xmx6g" SelectVariants -V ~{updated_input_vcf} \
        -L ~{vcf_intervals} \
        ~{"-L " + intersecting_intervals + " --interval-set-rule INTERSECTION"} \
        --exclude-filtered \
        --sites-only-vcf-output \
        --select-type-to-include SNP \
        --restrict-alleles-to BIALLELIC \
        --max-nocall-fraction 0.01 \
        -select "AF>0.001" \
        -O ~{output_filename}

        gatk IndexFeatureFile -I ~{output_filename} -O ~{output_filename_idx}
    >>>

    output {
        File vcf_sites_only="~{output_filename}"
        File vcf_sites_only_idx="~{output_filename_idx}"
    }

    runtime {
        docker:"us.gcr.io/broad-gatk/gatk:4.2.0.0"
        memory: "12 GB"
        cpu: "4"
        disks: "local-disk 100 HDD"
    }
}

# TODO: This is a dupe task from select_samples_for_training.wdl
task merge_vcf_bgzs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        String output_name
    }

    String output_vcf = basename(output_name) + ".vcf.bgz"
    String output_vcf_idx = basename(output_vcf) + ".tbi"

    parameter_meta {
        input_vcfs: {localization_optional: true}
        input_vcf_indices: {localization_optional: true}
    }

    command <<<
        set -e
        echo "Install Google CLI"
        apt-get update
        apt-get -y install curl python3
        curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-397.0.0-linux-x86_64.tar.gz
        ls -la
        tar -zxf google-cloud-cli-397.0.0-linux-x86_64.tar.gz
        ./google-cloud-sdk/install.sh
        ./google-cloud-sdk/bin/gcloud init --skip-diagnostics
        source /cromwell_root/google-cloud-sdk/completion.bash.inc
        source /cromwell_root/google-cloud-sdk/path.bash.inc

        echo "Doing manual localization..."
        mkdir -p vcf
        cp ~{write_lines(input_vcfs)} manifest_vcf.txt
        cp ~{write_lines(input_vcf_indices)} manifest_tbi.txt
        cd vcf

        # IMPORTANT:  This assumes that all input filenames are unique
        cat ../manifest_vcf.txt | gsutil -m cp -I .
        cat ../manifest_tbi.txt | gsutil -m cp -I .
        cd ..
        ls -t ${PWD}/vcf/*.bgz | sort --version-sort > local_vcfs.txt
        ls -t ${PWD}/vcf/*.tbi | sort --version-sort > local_tbis.txt

        head -2 local_vcfs.txt

        ## Tried to list the files directly on the command line, but this may have issues, so we have in a file of
        # filenames (fofn)
        echo "Concatenating..."
        date
        bcftools concat --naive -o tmp_~{output_vcf} -f local_vcfs.txt
        date
        echo "Starting sort... "
        date
        bcftools sort -o ~{output_vcf} tmp_~{output_vcf}
        date
        echo "Starting index... "
        bcftools index -t ~{output_vcf}
    >>>

    runtime {
        docker: "mgibio/bcftools-cwl:1.12"
        memory: "100 GB"
        cpu: "16"
        disks: "local-disk 1500 HDD"
        bootDiskSizeGb: 1500
    }

    output {
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
    }
}

# Output will be a bgz
task filter_by_sites_only {
    input {
        File vcf
        File vcf_idx
        File sites_only_vcf
        File sites_only_vcf_idx
        String output_name
        String id
        String? service_account_json
    }

    String output_filename = basename(output_name) + "." + id + ".vcf.bgz"
    String output_filename_idx = output_filename + ".tbi"
    String has_service_account_file = if (defined(service_account_json)) then 'true' else 'false'
    String service_account_basename_pre = if (defined(service_account_json)) then service_account_json else ''
    String service_account_basename = basename(service_account_basename_pre)
    String input_vcf_basename = basename(vcf)
    String updated_input_vcf = if (defined(service_account_json)) then input_vcf_basename else vcf

    parameter_meta {
        vcf: {localization_optional: true}
        vcf_idx: {localization_optional: true}
    }

    command <<<

        if [ ~{has_service_account_file} = 'true' ]; then
        gsutil cp ~{service_account_json} ~{service_account_basename}
        export GOOGLE_APPLICATION_CREDENTIALS=~{service_account_basename}
        gcloud auth activate-service-account --key-file='~{service_account_basename}'
        gsutil cp ~{vcf} .
        gsutil cp ~{vcf_idx} .
        fi

        gatk SelectVariants -V ~{updated_input_vcf} -L ~{sites_only_vcf}  -O ~{output_filename}
        gatk IndexFeatureFile -I ~{output_filename}
    >>>

    output {
        File filtered_vcf = "~{output_filename}"
        File filtered_vcf_idx = "~{output_filename_idx}"
    }

    runtime {
        docker:"us.gcr.io/broad-gatk/gatk:4.2.0.0"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
}

task intersect_vcfs_as_sites_only {
    input {
        File vcf1
        File vcf1_index
        File vcf2
        File vcf2_index
        String output1_basename
    }
    String output1_filename = output1_basename + ".vcf.bgz"
    String output1_filename_idx = output1_filename + ".tbi"

    command <<<
        gatk  SelectVariants \
        -V ~{vcf1} \
        -L ~{vcf2} \
        --sites-only-vcf-output \
        -O ~{output1_filename}

        gatk IndexFeatureFile -I ~{output1_filename} -O ~{output1_filename_idx}

        egrep -v "^\#" ~{output1_filename} | wc -l > variant_count_1.txt
    >>>
    output {
        File intersected_vcf = "~{output1_filename}"
        File intersected_vcf_idx = "~{output1_filename_idx}"
        File variant_count = "variant_count_1.txt"
    }

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.2.0.0"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 500 HDD"
    }
}