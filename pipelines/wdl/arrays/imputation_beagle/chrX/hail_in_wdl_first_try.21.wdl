version 1.0
#
# Given a VDS and a bed file, render a VCF (sharded by chromosome).
#   All bed files referenced in this WDL are UCSC bed files (as opposed to PLINK bed files).
#
# This has not been tested on any reference other than hg38.
# Inputs:
#
#         ## ANALYSIS PARAMETERS
#        #  ie, parameters that go to the Hail python code (submission_script below)
#        String vds_url
#
#        # Genomic region for the output VCFs to cover
#        File bed_url
#
#        # VCF Header that will be used in the output
#        File vcf_header_url
#
#        # Contigs of interest.  If a contig is present in the bed file, but not in this list, the contig will be ignored.
#        #   In other words, this is a contig level intersection with the bed file.
#        #     This list of contigs that must be present in the reference.  Each contig will be processed separately (shard)
#        # This list should be ordered.  Eg, ["chr21", "chr22"]
#        Array[String] contigs
#
#        # String used in construction of output filename
#        #  Cannot contain any special characters, ie, characters must be alphanumeric or "_"
#        String prefix
#
#        ## CLUSTER PARAMETERS
#        # Number of workers (per shard) to use in the Hail cluster.
#        Int num_workers
#
#        # The Google project ID information is necessary when spinning up dataproc.
#        #  ie, terra-<hex>
#        # This must match the workspace that this workflow is being run from.
#        #     eg, "terra-491d5f31"
#        String gcs_project
#
#        # Set to 'subnetwork' if running in Terra Cromwell
#        String gcs_subnetwork_name='subnetwork'
#
#        # The script that is run on the cluster
#        #  See filter_VDS_and_shard_by_contig.py for an example.
#        File submission_script
#
#        # Set to "us-central1" if running in Terra Cromwell
#        String region = "us-central1"
#
#        ## VM PARAMETERS
#        # Please note that there is a RuntimeAttr struct and a task parameter that can be used to override the defaults
#        #  of the VM.  These are task parameters.
#        #  However, since this can be a lightweight VM, overriding is unlikely to be necessary.
#
#        # The docker to be used on the VM.  This will need both Hail and Google Cloud SDK installed.
#        String hail_docker="us.gcr.io/broad-dsde-methods/lichtens/hail_dataproc_wdl:1.0"
#
# Important notes:
#   - Hail will save the VCFs in the cloud.  You will need to provide this storage space.  In other words, the runtime
#     parameters must have enough storage space to support a single contig
#   - This WDL script is still dependent on the python/Hail script that it calls.  You will see this when the parameters
#    are passed into the script.
#   - This WDL is boilerplate, except for input parameters, output parameters, and where marked in the main task.
#   - We HIGHLY recommend that the WDL do NOT run on a preemptible VM
#    (reminder, this is a single VM that spins up the dataproc cluster and submits jobs -- it is not doing any of the
#     actual computation.  In other words, it does not need to be a heavy machine.)
#     In other words, always set `preemptible_tries` to zero (default).
#
struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow RunAoUAnvilMergeFilterAndQc {
    ### Change here:  You will need to specify all parameters (both analysis and runtime) that need to go to the
    # cluster, VM spinning up the cluster, and the script being run on the cluster.
    input {

        ## ANALYSIS PARAMETERS
        #  ie, parameters that go to the Hail python code (submission_script below)

        # Input files/paths
        String aou_vds_url
        String anvil_vds_url
        String bed_url

        # Contigs of interest.  If a contig is present in the bed file, but not in this list, the contig will be ignored.
        #   In other words, this is a contig level intersection with the bed file.
        #     This list of contigs that must be present in the reference.  Each contig will be processed separately (shard)
        # This list should be ordered.  Eg, ["chr21", "chr22"]
        Array[String] contigs

        # String used in construction of output filename
        #  Cannot contain any special characters, ie, characters must be alphanumeric or "_"
        String prefix

        ## CLUSTER PARAMETERS

        # The Google project ID information is necessary when spinning up dataproc.
        #  ie, terra-<hex>
        # This must match the workspace that this workflow is being run from.
        #     eg, "terra-491d5f31"
        String gcs_project

        # Set to 'subnetwork' if running in Terra Cromwell
        String gcs_subnetwork_name = 'subnetwork'

        # The script that is run on the cluster
        File submission_script_step1
        File submission_script_step2
        File submission_script_step3
        File submission_script_step3b
        File submission_script_step4
        File submission_script_step5a
        File submission_script_step5b

        # Bucket used to stage input and output files, likely the workspace bucket. 
        #  Include "gs://" prefix and any subdirectory structure, but no trailing "/"
        String output_bucket_path

        # Set to "us-central1" if running in Terra Cromwell
        String region = "us-central1"

        ## VM PARAMETERS
        # Please note that there is a RuntimeAttr struct and a task parameter that can be used to override the defaults
        #  of the VM.  These are task parameters.
        #  However, since this can be a lightweight VM, overriding is unlikely to be necessary.

        # The docker to be used on the VM.  This will need both Hail and Google Cloud SDK installed.
        String hail_docker = "gcr.io/broad-dsde-methods/aou-auxiliary/hail_dataproc_wdl:0.2.130"
    }

    # because people can't read directions
    String output_bucket_path_with_trailing_slash = sub(output_bucket_path, "/$", "") + "/"

    scatter (contig in contigs) {

        call Step1 {
            input:
                input_aou_vds_url = aou_vds_url,
                contig = contig,
                bed_url = bed_url,
                prefix = prefix,
                gcs_project = gcs_project,
                gcs_subnetwork_name = gcs_subnetwork_name,
                submission_script = submission_script_step1,
                output_bucket = output_bucket_path_with_trailing_slash,
                hail_docker = hail_docker,
                region = region
        }

        call Step2 {
            input: 
                input_anvil_vds_url = anvil_vds_url,
                contig = contig,
                bed_url = bed_url,
                prefix = prefix,
                gcs_project = gcs_project,
                gcs_subnetwork_name = gcs_subnetwork_name,
                submission_script = submission_script_step2,
                output_bucket = output_bucket_path_with_trailing_slash,
                hail_docker = hail_docker,
                region = region
        }

        call Step3 {
            input: 
                input_aou_vds_url = aou_vds_url,
                input_anvil_vds_url = anvil_vds_url,
                input_aou_vcf_url = Step1.aou_vcf,
                input_anvil_mt_url = Step2.anvil_mt,
                contig = contig,
                bed_url = bed_url,
                prefix = prefix,
                gcs_project = gcs_project,
                gcs_subnetwork_name = gcs_subnetwork_name,
                submission_script = submission_script_step3,
                output_bucket = output_bucket_path_with_trailing_slash,
                hail_docker = hail_docker,
                region = region
        }

        call Step3b {
            input:
                input_anvil_filtered_vds_url = Step3.anvil_filtered_vds,
                contig = contig,
                bed_url = bed_url,
                prefix = prefix,
                gcs_project = gcs_project,
                gcs_subnetwork_name = gcs_subnetwork_name,
                submission_script = submission_script_step3b,
                output_bucket = output_bucket_path_with_trailing_slash,
                hail_docker = hail_docker,
                region = region
        }

        call Step4 {
            input: 
                input_aou_filtered_vds_url = Step3.aou_filtered_vds,
                input_anvil_filtered_vds_url = Step3b.anvil_renamed_vds,
                contig = contig,
                bed_url = bed_url,
                prefix = prefix,
                gcs_project = gcs_project,
                gcs_subnetwork_name = gcs_subnetwork_name,
                submission_script = submission_script_step4,
                output_bucket = output_bucket_path_with_trailing_slash,
                hail_docker = hail_docker,
                region = region
        }

        call Step5a {
            input: 
                input_aou_vcf_url = Step1.aou_vcf,
                input_anvil_mt_url = Step2.anvil_mt,
                input_combined_vds_url = Step4.combined_vds,
                contig = contig,
                bed_url = bed_url,
                prefix = prefix,
                gcs_project = gcs_project,
                gcs_subnetwork_name = gcs_subnetwork_name,
                submission_script = submission_script_step5a,
                output_bucket = output_bucket_path_with_trailing_slash,
                hail_docker = hail_docker,
                region = region
        }

        call Step5b {
            input:
                input_filtered_vds_url = Step5a.filtered_vds,
                contig = contig,
                bed_url = bed_url,
                prefix = prefix,
                gcs_project = gcs_project,
                gcs_subnetwork_name = gcs_subnetwork_name,
                submission_script = submission_script_step5b,
                output_bucket = output_bucket_path_with_trailing_slash,
                hail_docker = hail_docker,
                region = region
        }
    }

    output {
        # Note we can remove any of these if they're purely intermediates and have no value as final outputs

        # Step 1
        Array[String] aou_vcfs = Step1.aou_vcf
        Array[String] aou_vcf_headers = Step1.aou_vcf_header
        Array[String] step1_reports = Step1.report

        # Step 2
        Array[String] anvil_mts = Step2.anvil_mt
        Array[String] step2_reports = Step2.report

        # Step 3
        Array[String] aou_filtered_vdss = Step3.aou_filtered_vds
        Array[String] anvil_filtered_vdss = Step3.anvil_filtered_vds
        Array[String] step3_reports = Step3.report

        # Step 4
        Array[String] combined_vdss = Step4.combined_vds
        Array[String] step4_reports = Step4.report

        # Step 5a
        Array[String] filtered_vds = Step5a.filtered_vds
        Array[String] step5a_reports = Step5a.report

        # Step 5b
        Array[String] final_vcf = Step5b.final_vcf
        Array[String] step5b_reports = Step5b.report

    }
}

task Step1 {
    input {
        # You must treat a VDS as a String, since it is a directory and not a single file
        String input_aou_vds_url
        String bed_url

        # Cannot make this localization_optional
        File submission_script
        String output_bucket

        # contig must be in the reference
        String contig
        String prefix

        # dataproc params
        String gcs_project
        String region = "us-central1"
        String master_machine_type = "n1-highmem-32"
        Float master_memory_fraction = 0.8
        String worker_machine_type = "n1-highmem-4"
        Int num_workers = 2
        Int num_preemptible_workers = 50
        Int time_to_live_minutes = 2880 # two days
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name

        String hail_docker
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 6.5,
                                      disk_gb: 150,
                                      cpu_cores: 1,
                                      preemptible_tries: 0,
                                      max_retries: 0,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    # define output file urls with file names
    String output_aou_vcf_url        = output_bucket + contig + ".step1.aou.vcf.bgz"
    String output_aou_vcf_header_url = output_bucket + contig + ".step1.aou.vcf.header"
    String output_report_url         = output_bucket + contig + ".step1.report"

    command <<<
        set -euxo pipefail

        gcloud config list account --format "value(core.account)" 1> account.txt

        #### TEST:  Make sure that this docker image is configured for python3
        if which python3 > /dev/null 2>&1; then
            pt3="$(which python3)"
            echo "** python3 located at $pt3"
            echo "** magic: $(file $pt3)"
            echo "** Version info:"
            echo "$(python3 -V)"
            echo "** -c test"
            python3 -c "print('hello world')"
        else
            echo "!! No 'python3' in path."
            exit 1
        fi
        #### END TEST

        python3 <<EOF
        print("Running python code...")
        import hail as hl
        import os
        import uuid
        from google.cloud import dataproc_v1 as dataproc

        # Must match pattern (?:[a-z](?:[-a-z0-9]{0,49}[a-z0-9])?)
        cluster_name = f'~{prefix}-~{contig}-hail-step1-{str(uuid.uuid4())[0:13]}'.lower()

        # Must be local filepath
        script_path = "~{submission_script}"

        with open("account.txt", "r") as account_file:
            account = account_file.readline().strip()
        print("account: " + account)

        try:
            cluster_start_cmd = "hailctl dataproc start --master-machine-type {} --master-memory-fraction ~{master_memory_fraction} --worker-machine-type {} --num-workers ~{num_workers} --num-preemptible-workers ~{num_preemptible_workers} --region {} --project {} --service-account {} --num-master-local-ssds 1 --num-worker-local-ssds 1 --max-idle=60m --max-age=~{time_to_live_minutes}m --subnet={} {}".format("~{master_machine_type}", "~{worker_machine_type}", "~{region}", "~{gcs_project}", account, "projects/~{gcs_project}/regions/~{region}/subnetworks/~{gcs_subnetwork_name}", cluster_name)
            print("Starting cluster...")
            print(cluster_start_cmd)
            f = os.popen(cluster_start_cmd)
            f.read()
            if (f.close() != None):
                raise Exception("Failed to start cluster sucessfully")

            cluster_client = dataproc.ClusterControllerClient(
                client_options={"api_endpoint": f"~{region}-dataproc.googleapis.com:443"}
            )

            for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
                if cluster.cluster_name == cluster_name:
                    cluster_temp_bucket = cluster.config.temp_bucket

                    #### THIS IS WHERE YOU CALL YOUR SCRIPT AND COPY THE OUTPUT LOCALLY (so that it can get back into WDL-space)
                    submit_cmd = f'''gcloud dataproc jobs submit pyspark {script_path} \
                    --cluster={cluster_name} --project ~{gcs_project} --region=~{region} --account {account} --driver-log-levels root=WARN -- \
                    --input_aou_vds_url ~{input_aou_vds_url} \
                    --output_aou_vcf_url ~{output_aou_vcf_url} \
                    --output_aou_vcf_header_url ~{output_aou_vcf_header_url} \
                    --output_report_url ~{output_report_url} \
                    --contig ~{contig} \
                    --bed_url ~{bed_url} \
                    --temp_bucket gs://{cluster_temp_bucket}/{cluster_name}'''

                    print("Running: " + submit_cmd)
                    f = os.popen(submit_cmd)
                    f.read()
                    if (f.close() != None):
                        raise Exception("Failed to submit cluster job sucessfully")
                    ###########

                    break

        except Exception as e:
            print(e)
            raise
        finally:
            print(f'Stopping cluster: {cluster_name}')
            os.popen("gcloud dataproc clusters delete --project {} --region {} --account {} {}".format("~{gcs_project}", "~{region}", account, cluster_name)).read()

        EOF

        echo "Complete"
    >>>

    output {
        String aou_vcf = output_aou_vcf_url
        String aou_vcf_header = output_aou_vcf_header_url
        String report = output_report_url
    }

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}

task Step2 {
    input {
        # You must treat a VDS as a String, since it is a directory and not a single file
        String input_anvil_vds_url
        String bed_url

        # Cannot make this localization_optional
        File submission_script
        String output_bucket

        # contig must be in the reference
        String contig
        String prefix

        # dataproc params
        String gcs_project
        String region = "us-central1"
        String master_machine_type = "n1-highmem-32"
        Float master_memory_fraction = 0.8
        String worker_machine_type = "n1-highmem-4"
        Int num_workers = 2
        Int num_preemptible_workers = 50
        Int time_to_live_minutes = 2880 # two days
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name

        String hail_docker
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 6.5,
                                      disk_gb: 150,
                                      cpu_cores: 1,
                                      preemptible_tries: 0,
                                      max_retries: 0,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    # define output file urls with file names
    String output_anvil_mt_url = output_bucket + contig + ".step2.anvil.mt"
    String output_report_url   = output_bucket + contig + ".step2.report"

    command <<<
        set -euxo pipefail

        gcloud config list account --format "value(core.account)" 1> account.txt

        #### TEST:  Make sure that this docker image is configured for python3
        if which python3 > /dev/null 2>&1; then
            pt3="$(which python3)"
            echo "** python3 located at $pt3"
            echo "** magic: $(file $pt3)"
            echo "** Version info:"
            echo "$(python3 -V)"
            echo "** -c test"
            python3 -c "print('hello world')"
        else
            echo "!! No 'python3' in path."
            exit 1
        fi
        #### END TEST

        python3 <<EOF
        print("Running python code...")
        import hail as hl
        import os
        import uuid
        from google.cloud import dataproc_v1 as dataproc

        # Must match pattern (?:[a-z](?:[-a-z0-9]{0,49}[a-z0-9])?)
        cluster_name = f'~{prefix}-~{contig}-hail-step2-{str(uuid.uuid4())[0:13]}'.lower()

        # Must be local filepath
        script_path = "~{submission_script}"

        with open("account.txt", "r") as account_file:
            account = account_file.readline().strip()
        print("account: " + account)

        try:
            cluster_start_cmd = "hailctl dataproc start --master-machine-type {} --master-memory-fraction ~{master_memory_fraction} --worker-machine-type {} --num-workers ~{num_workers} --num-preemptible-workers ~{num_preemptible_workers} --region {} --project {} --service-account {} --num-master-local-ssds 1 --num-worker-local-ssds 1 --max-idle=60m --max-age=~{time_to_live_minutes}m --requester-pays-allow-all --subnet={} {}".format("~{master_machine_type}", "~{worker_machine_type}", "~{region}", "~{gcs_project}", account, "projects/~{gcs_project}/regions/~{region}/subnetworks/~{gcs_subnetwork_name}", cluster_name)
            print("Starting cluster...")
            print(cluster_start_cmd)
            f = os.popen(cluster_start_cmd)
            f.read()
            if (f.close() != None):
                raise Exception("Failed to start cluster sucessfully")

            cluster_client = dataproc.ClusterControllerClient(
                client_options={"api_endpoint": f"~{region}-dataproc.googleapis.com:443"}
            )

            for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
                if cluster.cluster_name == cluster_name:
                    cluster_temp_bucket = cluster.config.temp_bucket

                    #### THIS IS WHERE YOU CALL YOUR SCRIPT AND COPY THE OUTPUT LOCALLY (so that it can get back into WDL-space)
                    submit_cmd = f'''gcloud dataproc jobs submit pyspark {script_path} \
                    --cluster={cluster_name} --project ~{gcs_project} --region=~{region} --account {account} --driver-log-levels root=WARN -- \
                    --input_anvil_vds_url ~{input_anvil_vds_url} \
                    --output_anvil_mt_url ~{output_anvil_mt_url} \
                    --output_report_url ~{output_report_url} \
                    --contig ~{contig} \
                    --bed_url ~{bed_url} \
                    --temp_bucket gs://{cluster_temp_bucket}/{cluster_name}'''
                    print("Running: " + submit_cmd)
                    f = os.popen(submit_cmd)
                    f.read()
                    if (f.close() != None):
                        raise Exception("Failed to submit cluster job sucessfully")
                    ###########

                    break

        except Exception as e:
            print(e)
            raise
        finally:
            print(f'Stopping cluster: {cluster_name}')
            os.popen("gcloud dataproc clusters delete --project {} --region {} --account {} {}".format("~{gcs_project}", "~{region}", account, cluster_name)).read()

        EOF

        echo "Complete"
    >>>

    output {
        String anvil_mt = output_anvil_mt_url
        String report = output_report_url
    }

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}

task Step3 {
    input {
        # You must treat a VDS as a String, since it is a directory and not a single file
        String input_aou_vds_url
        String input_anvil_vds_url

        String input_aou_vcf_url     # output from Step 1
        String input_anvil_mt_url    # output from Step 2
        String bed_url

        # Cannot make this localization_optional
        File submission_script
        String output_bucket

        # contig must be in the reference
        String contig
        String prefix

        # dataproc params
        String gcs_project
        String region = "us-central1"
        String master_machine_type = "n1-highmem-32"
        Float master_memory_fraction = 0.8
        String worker_machine_type = "n1-highmem-4"
        Int num_workers = 2
        Int num_preemptible_workers = 50
        Int time_to_live_minutes = 2880 # two days
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name

        String hail_docker
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 6.5,
                                      disk_gb: 150,
                                      cpu_cores: 1,
                                      preemptible_tries: 0,
                                      max_retries: 0,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    # define output file urls with file names
    String output_aou_filtered_vds_url   = output_bucket + contig + ".step3.aou.filtered.vds"
    String output_anvil_filtered_vds_url = output_bucket + contig + ".step3.anvil.filtered.vds"
    String output_report_url             = output_bucket + contig + ".step3.report"

    command <<<
        set -euxo pipefail

        gcloud config list account --format "value(core.account)" 1> account.txt

        #### TEST:  Make sure that this docker image is configured for python3
        if which python3 > /dev/null 2>&1; then
            pt3="$(which python3)"
            echo "** python3 located at $pt3"
            echo "** magic: $(file $pt3)"
            echo "** Version info:"
            echo "$(python3 -V)"
            echo "** -c test"
            python3 -c "print('hello world')"
        else
            echo "!! No 'python3' in path."
            exit 1
        fi
        #### END TEST

        python3 <<EOF
        print("Running python code...")
        import hail as hl
        import os
        import uuid
        from google.cloud import dataproc_v1 as dataproc

        # Must match pattern (?:[a-z](?:[-a-z0-9]{0,49}[a-z0-9])?)
        cluster_name = f'~{prefix}-~{contig}-hail-step3-{str(uuid.uuid4())[0:13]}'.lower()

        # Must be local filepath
        script_path = "~{submission_script}"

        with open("account.txt", "r") as account_file:
            account = account_file.readline().strip()
        print("account: " + account)

        try:
            cluster_start_cmd = "hailctl dataproc start --master-machine-type {} --master-memory-fraction ~{master_memory_fraction} --worker-machine-type {} --num-workers ~{num_workers} --num-preemptible-workers ~{num_preemptible_workers} --region {} --project {} --service-account {} --num-master-local-ssds 1 --num-worker-local-ssds 1 --max-idle=60m --max-age=~{time_to_live_minutes}m --requester-pays-allow-all --subnet={} {}".format("~{master_machine_type}", "~{worker_machine_type}", "~{region}", "~{gcs_project}", account, "projects/~{gcs_project}/regions/~{region}/subnetworks/~{gcs_subnetwork_name}", cluster_name)
            print("Starting cluster...")
            print(cluster_start_cmd)
            f = os.popen(cluster_start_cmd)
            f.read()
            if (f.close() != None):
                raise Exception("Failed to start cluster sucessfully")

            cluster_client = dataproc.ClusterControllerClient(
                client_options={"api_endpoint": f"~{region}-dataproc.googleapis.com:443"}
            )

            for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
                if cluster.cluster_name == cluster_name:
                    cluster_temp_bucket = cluster.config.temp_bucket

                    #### THIS IS WHERE YOU CALL YOUR SCRIPT AND COPY THE OUTPUT LOCALLY (so that it can get back into WDL-space)
                    submit_cmd = f'''gcloud dataproc jobs submit pyspark {script_path} \
                    --cluster={cluster_name} --project ~{gcs_project} --region=~{region} --account {account} --driver-log-levels root=WARN -- \
                    --input_aou_vcf_url ~{input_aou_vcf_url} \
                    --input_anvil_mt_url ~{input_anvil_mt_url} \
                    --input_aou_vds_url ~{input_aou_vds_url} \
                    --input_anvil_vds_url ~{input_anvil_vds_url} \
                    --output_aou_filtered_vds_url ~{output_aou_filtered_vds_url} \
                    --output_anvil_filtered_vds_url ~{output_anvil_filtered_vds_url} \
                    --output_report_url ~{output_report_url} \
                    --contig ~{contig} \
                    --bed_url ~{bed_url} \
                    --temp_bucket gs://{cluster_temp_bucket}/{cluster_name}'''

                    print("Running: " + submit_cmd)
                    f = os.popen(submit_cmd)
                    f.read()
                    if (f.close() != None):
                        raise Exception("Failed to submit cluster job sucessfully")
                    ###########

                    break

        except Exception as e:
            print(e)
            raise
        finally:
            print(f'Stopping cluster: {cluster_name}')
            os.popen("gcloud dataproc clusters delete --project {} --region {} --account {} {}".format("~{gcs_project}", "~{region}", account, cluster_name)).read()

        EOF

        echo "Complete"
    >>>

    output {
        String aou_filtered_vds = output_aou_filtered_vds_url
        String anvil_filtered_vds = output_anvil_filtered_vds_url
        String report = output_report_url
    }

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}

task Step3b {
    input {
        # You must treat a VDS as a String, since it is a directory and not a single file
        String input_anvil_filtered_vds_url # output from Step 3
        String bed_url

        # Cannot make this localization_optional
        File submission_script
        String output_bucket

        # contig must be in the reference
        String contig
        String prefix

        # dataproc params
        String gcs_project
        String region = "us-central1"
        String master_machine_type = "n1-highmem-32"
        Float master_memory_fraction = 0.8
        String worker_machine_type = "n1-highmem-4"
        Int num_workers = 2
        Int num_preemptible_workers = 50
        Int time_to_live_minutes = 2880 # two days
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name

        String hail_docker
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 6.5,
                                      disk_gb: 150,
                                      cpu_cores: 1,
                                      preemptible_tries: 0,
                                      max_retries: 0,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    # define output file urls with file names
    String output_anvil_renamed_vds_url = output_bucket + contig + ".step3b.anvil.renamed.vds"
    String output_report_url            = output_bucket + contig + ".step3b.report"

    command <<<
        set -euxo pipefail

        gcloud config list account --format "value(core.account)" 1> account.txt

        #### TEST:  Make sure that this docker image is configured for python3
        if which python3 > /dev/null 2>&1; then
            pt3="$(which python3)"
            echo "** python3 located at $pt3"
            echo "** magic: $(file $pt3)"
            echo "** Version info:"
            echo "$(python3 -V)"
            echo "** -c test"
            python3 -c "print('hello world')"
        else
            echo "!! No 'python3' in path."
            exit 1
        fi
        #### END TEST

        python3 <<EOF
        print("Running python code...")
        import hail as hl
        import os
        import uuid
        from google.cloud import dataproc_v1 as dataproc

        # Must match pattern (?:[a-z](?:[-a-z0-9]{0,49}[a-z0-9])?)
        cluster_name = f'~{prefix}-~{contig}-hail-step3b-{str(uuid.uuid4())[0:13]}'.lower()

        # Must be local filepath
        script_path = "~{submission_script}"

        with open("account.txt", "r") as account_file:
            account = account_file.readline().strip()
        print("account: " + account)

        try:
            cluster_start_cmd = "hailctl dataproc start --master-machine-type {} --master-memory-fraction ~{master_memory_fraction} --worker-machine-type {} --num-workers ~{num_workers} --num-preemptible-workers ~{num_preemptible_workers} --region {} --project {} --service-account {} --num-master-local-ssds 1 --num-worker-local-ssds 1 --max-idle=60m --max-age=~{time_to_live_minutes}m --requester-pays-allow-all --subnet={} {}".format("~{master_machine_type}", "~{worker_machine_type}", "~{region}", "~{gcs_project}", account, "projects/~{gcs_project}/regions/~{region}/subnetworks/~{gcs_subnetwork_name}", cluster_name)
            print("Starting cluster...")
            print(cluster_start_cmd)
            f = os.popen(cluster_start_cmd)
            f.read()
            if (f.close() != None):
                raise Exception("Failed to start cluster sucessfully")

            cluster_client = dataproc.ClusterControllerClient(
                client_options={"api_endpoint": f"~{region}-dataproc.googleapis.com:443"}
            )

            for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
                if cluster.cluster_name == cluster_name:
                    cluster_temp_bucket = cluster.config.temp_bucket

                    #### THIS IS WHERE YOU CALL YOUR SCRIPT AND COPY THE OUTPUT LOCALLY (so that it can get back into WDL-space)
                    submit_cmd = f'''gcloud dataproc jobs submit pyspark {script_path} \
                    --cluster={cluster_name} --project ~{gcs_project} --region=~{region} --account {account} --driver-log-levels root=WARN -- \
                    --input_anvil_filtered_vds_url ~{input_anvil_filtered_vds_url} \
                    --output_anvil_renamed_vds_url ~{output_anvil_renamed_vds_url} \
                    --output_report_url ~{output_report_url} \
                    --contig ~{contig} \
                    --bed_url ~{bed_url} \
                    --temp_bucket gs://{cluster_temp_bucket}/{cluster_name}'''

                    print("Running: " + submit_cmd)
                    f = os.popen(submit_cmd)
                    f.read()
                    if (f.close() != None):
                        raise Exception("Failed to submit cluster job sucessfully")
                    ###########

                    break

        except Exception as e:
            print(e)
            raise
        finally:
            print(f'Stopping cluster: {cluster_name}')
            os.popen("gcloud dataproc clusters delete --project {} --region {} --account {} {}".format("~{gcs_project}", "~{region}", account, cluster_name)).read()

        EOF

        echo "Complete"
    >>>

    output {
        String anvil_renamed_vds = output_anvil_renamed_vds_url
        String report = output_report_url
    }

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}

task Step4 {
    input {
        # You must treat a VDS as a String, since it is a directory and not a single file
        String input_aou_filtered_vds_url       # output from Step 3
        String input_anvil_filtered_vds_url     # output from Step 3
        String bed_url

        # Cannot make this localization_optional
        File submission_script
        String output_bucket

        # contig must be in the reference
        String contig
        String prefix

        # dataproc params
        String gcs_project
        String region = "us-central1"
        String master_machine_type = "n1-highmem-32"
        Float master_memory_fraction = 0.8
        String worker_machine_type = "n1-highmem-4"
        Int num_workers = 2
        Int num_preemptible_workers = 50
        Int time_to_live_minutes = 10080 # seven days
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name

        String hail_docker
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 6.5,
                                      disk_gb: 150,
                                      cpu_cores: 1,
                                      preemptible_tries: 0,
                                      max_retries: 0,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    # define output file urls with file names
    String output_combined_vds_url = output_bucket + contig + ".step4.aou-anvil.combined.vds"
    String output_report_url       = output_bucket + contig + ".step4.report"

    command <<<
        set -euxo pipefail

        gcloud config list account --format "value(core.account)" 1> account.txt

        #### TEST:  Make sure that this docker image is configured for python3
        if which python3 > /dev/null 2>&1; then
            pt3="$(which python3)"
            echo "** python3 located at $pt3"
            echo "** magic: $(file $pt3)"
            echo "** Version info:"
            echo "$(python3 -V)"
            echo "** -c test"
            python3 -c "print('hello world')"
        else
            echo "!! No 'python3' in path."
            exit 1
        fi
        #### END TEST

        python3 <<EOF
        print("Running python code...")
        import hail as hl
        import os
        import uuid
        from google.cloud import dataproc_v1 as dataproc

        # Must match pattern (?:[a-z](?:[-a-z0-9]{0,49}[a-z0-9])?)
        cluster_name = f'~{prefix}-~{contig}-hail-step4-{str(uuid.uuid4())[0:13]}'.lower()

        # Must be local filepath
        script_path = "~{submission_script}"

        with open("account.txt", "r") as account_file:
            account = account_file.readline().strip()
        print("account: " + account)

        try:
            cluster_start_cmd = "hailctl dataproc start --master-machine-type {} --master-memory-fraction ~{master_memory_fraction} --worker-machine-type {} --num-workers ~{num_workers} --num-preemptible-workers ~{num_preemptible_workers} --region {} --project {} --service-account {} --num-master-local-ssds 1 --num-worker-local-ssds 1 --max-idle=60m --max-age=~{time_to_live_minutes}m --subnet={} {}".format("~{master_machine_type}", "~{worker_machine_type}", "~{region}", "~{gcs_project}", account, "projects/~{gcs_project}/regions/~{region}/subnetworks/~{gcs_subnetwork_name}", cluster_name)
            print("Starting cluster...")
            print(cluster_start_cmd)
            f = os.popen(cluster_start_cmd)
            f.read()
            if (f.close() != None):
                raise Exception("Failed to start cluster sucessfully")

            cluster_client = dataproc.ClusterControllerClient(
                client_options={"api_endpoint": f"~{region}-dataproc.googleapis.com:443"}
            )

            for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
                if cluster.cluster_name == cluster_name:
                    cluster_temp_bucket = cluster.config.temp_bucket

                    #### THIS IS WHERE YOU CALL YOUR SCRIPT AND COPY THE OUTPUT LOCALLY (so that it can get back into WDL-space)
                    submit_cmd = f'''gcloud dataproc jobs submit pyspark {script_path} \
                    --cluster={cluster_name} --project ~{gcs_project} --region=~{region} --account {account} --driver-log-levels root=WARN -- \
                    --input_aou_filtered_vds_url ~{input_aou_filtered_vds_url} \
                    --input_anvil_filtered_vds_url ~{input_anvil_filtered_vds_url} \
                    --output_combined_vds_url ~{output_combined_vds_url} \
                    --output_report_url ~{output_report_url} \
                    --contig ~{contig} \
                    --bed_url ~{bed_url} \
                    --temp_bucket gs://{cluster_temp_bucket}/{cluster_name}'''

                    print("Running: " + submit_cmd)
                    f = os.popen(submit_cmd)
                    f.read()
                    if (f.close() != None):
                        raise Exception("Failed to submit cluster job sucessfully")
                    ###########

                    break

        except Exception as e:
            print(e)
            raise
        finally:
            print(f'Stopping cluster: {cluster_name}')
            os.popen("gcloud dataproc clusters delete --project {} --region {} --account {} {}".format("~{gcs_project}", "~{region}", account, cluster_name)).read()

        EOF

        echo "Complete"
    >>>

    output {
        String combined_vds = output_combined_vds_url
        String report = output_report_url
    }

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}

task Step5a {
    input {
        # You must treat a VDS as a String, since it is a directory and not a single file
        String input_aou_vcf_url            # output from Step 1
        String input_anvil_mt_url           # output from Step 2
        String input_combined_vds_url       # output from Step 4
        String bed_url

        # Cannot make this localization_optional
        File submission_script
        String output_bucket

        # contig must be in the reference
        String contig
        String prefix

        # dataproc params
        String gcs_project
        String region = "us-central1"
        String master_machine_type = "n1-highmem-64"
        Float master_memory_fraction = 0.8
        String worker_machine_type = "n1-highmem-8"
        Int num_workers = 2
        Int num_preemptible_workers = 50
        Int time_to_live_minutes = 2880 # two days
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name

        String hail_docker
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 6.5,
                                      disk_gb: 150,
                                      cpu_cores: 1,
                                      preemptible_tries: 0,
                                      max_retries: 0,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    # define output file urls with file names
    String output_filtered_vds_url        = output_bucket + "aou-anvil.filtered.step5a." + contig + ".filtered.vds"
    String output_report_url              = output_bucket + contig + ".step5a.report"

    command <<<
        set -euxo pipefail

        gcloud config list account --format "value(core.account)" 1> account.txt

        #### TEST:  Make sure that this docker image is configured for python3
        if which python3 > /dev/null 2>&1; then
            pt3="$(which python3)"
            echo "** python3 located at $pt3"
            echo "** magic: $(file $pt3)"
            echo "** Version info:"
            echo "$(python3 -V)"
            echo "** -c test"
            python3 -c "print('hello world')"
        else
            echo "!! No 'python3' in path."
            exit 1
        fi
        #### END TEST

        python3 <<EOF
        print("Running python code...")
        import hail as hl
        import os
        import uuid
        from google.cloud import dataproc_v1 as dataproc

        # Must match pattern (?:[a-z](?:[-a-z0-9]{0,49}[a-z0-9])?)
        cluster_name = f'~{prefix}-~{contig}-hail-step5a-{str(uuid.uuid4())[0:13]}'.lower()

        # Must be local filepath
        script_path = "~{submission_script}"

        with open("account.txt", "r") as account_file:
            account = account_file.readline().strip()
        print("account: " + account)

        try:
            cluster_start_cmd = "hailctl dataproc start --master-machine-type {} --master-memory-fraction ~{master_memory_fraction} --worker-machine-type {} --num-workers ~{num_workers} --num-preemptible-workers ~{num_preemptible_workers} --region {} --project {} --service-account {} --num-master-local-ssds 1 --num-worker-local-ssds 1 --max-idle=60m --max-age=~{time_to_live_minutes}m --subnet={} {}".format("~{master_machine_type}", "~{worker_machine_type}", "~{region}", "~{gcs_project}", account, "projects/~{gcs_project}/regions/~{region}/subnetworks/~{gcs_subnetwork_name}", cluster_name)
            print("Starting cluster...")
            print(cluster_start_cmd)
            f = os.popen(cluster_start_cmd)
            f.read()
            if (f.close() != None):
                raise Exception("Failed to start cluster sucessfully")

            cluster_client = dataproc.ClusterControllerClient(
                client_options={"api_endpoint": f"~{region}-dataproc.googleapis.com:443"}
            )

            for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
                if cluster.cluster_name == cluster_name:
                    cluster_temp_bucket = cluster.config.temp_bucket

                    #### THIS IS WHERE YOU CALL YOUR SCRIPT AND COPY THE OUTPUT LOCALLY (so that it can get back into WDL-space)
                    submit_cmd = f'''gcloud dataproc jobs submit pyspark {script_path} \
                    --cluster={cluster_name} --project ~{gcs_project} --region=~{region} --account {account} --driver-log-levels root=WARN -- \
                    --input_combined_vds_url ~{input_combined_vds_url} \
                    --input_anvil_mt_url ~{input_anvil_mt_url} \
                    --input_aou_vcf_url ~{input_aou_vcf_url} \
                    --output_filtered_vds_url ~{output_filtered_vds_url} \
                    --output_report_url ~{output_report_url} \
                    --contig ~{contig} \
                    --bed_url ~{bed_url} \
                    --temp_bucket gs://{cluster_temp_bucket}/{cluster_name}'''

                    print("Running: " + submit_cmd)
                    f = os.popen(submit_cmd)
                    f.read()
                    if (f.close() != None):
                        raise Exception("Failed to submit cluster job sucessfully")
                    ###########

                    break

        except Exception as e:
            print(e)
            raise
        finally:
            print(f'Stopping cluster: {cluster_name}')
            os.popen("gcloud dataproc clusters delete --project {} --region {} --account {} {}".format("~{gcs_project}", "~{region}", account, cluster_name)).read()

        EOF

        echo "Complete"
    >>>

    output {
        String filtered_vds = output_filtered_vds_url
        String report = output_report_url
    }

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}

task Step5b {
    input {
        # You must treat a VDS as a String, since it is a directory and not a single file
        String input_filtered_vds_url      # output from Step 5a
        String bed_url

        # Cannot make this localization_optional
        File submission_script
        String output_bucket

        # contig must be in the reference
        String contig
        String prefix

        # dataproc params
        String gcs_project
        String region = "us-central1"
        String master_machine_type = "n1-highmem-64"
        Float master_memory_fraction = 0.8
        String worker_machine_type = "n1-highmem-8"
        Int num_workers = 2
        Int num_preemptible_workers = 50
        Int time_to_live_minutes = 2880 # two days
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name

        String hail_docker
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 6.5,
                                      disk_gb: 150,
                                      cpu_cores: 1,
                                      preemptible_tries: 0,
                                      max_retries: 0,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    # define output file urls with file names
    String output_final_vcf_url           = output_bucket + "aou-anvil.final.step5b." + contig + ".vcf.bgz" # imputation wdl expects contig right before prefix
    String output_final_vcf_header_url    = output_bucket + "aou-anvil.final." + contig + ".step5b.vcf.header"
    String output_report_url              = output_bucket + contig + ".step5b.report"

    command <<<
        set -euxo pipefail

        gcloud config list account --format "value(core.account)" 1> account.txt

        #### TEST:  Make sure that this docker image is configured for python3
        if which python3 > /dev/null 2>&1; then
            pt3="$(which python3)"
            echo "** python3 located at $pt3"
            echo "** magic: $(file $pt3)"
            echo "** Version info:"
            echo "$(python3 -V)"
            echo "** -c test"
            python3 -c "print('hello world')"
        else
            echo "!! No 'python3' in path."
            exit 1
        fi
        #### END TEST

        python3 <<EOF
        print("Running python code...")
        import hail as hl
        import os
        import uuid
        from google.cloud import dataproc_v1 as dataproc

        # Must match pattern (?:[a-z](?:[-a-z0-9]{0,49}[a-z0-9])?)
        cluster_name = f'~{prefix}-~{contig}-hail-step5b-{str(uuid.uuid4())[0:13]}'.lower()

        # Must be local filepath
        script_path = "~{submission_script}"

        with open("account.txt", "r") as account_file:
            account = account_file.readline().strip()
        print("account: " + account)

        try:
            cluster_start_cmd = "hailctl dataproc start --master-machine-type {} --master-memory-fraction ~{master_memory_fraction} --worker-machine-type {} --num-workers ~{num_workers} --num-preemptible-workers ~{num_preemptible_workers} --region {} --project {} --service-account {} --num-master-local-ssds 1 --num-worker-local-ssds 1 --max-idle=60m --max-age=~{time_to_live_minutes}m --subnet={} {}".format("~{master_machine_type}", "~{worker_machine_type}", "~{region}", "~{gcs_project}", account, "projects/~{gcs_project}/regions/~{region}/subnetworks/~{gcs_subnetwork_name}", cluster_name)
            print("Starting cluster...")
            print(cluster_start_cmd)
            f = os.popen(cluster_start_cmd)
            f.read()
            if (f.close() != None):
                raise Exception("Failed to start cluster sucessfully")

            cluster_client = dataproc.ClusterControllerClient(
                client_options={"api_endpoint": f"~{region}-dataproc.googleapis.com:443"}
            )

            for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
                if cluster.cluster_name == cluster_name:
                    cluster_temp_bucket = cluster.config.temp_bucket

                    #### THIS IS WHERE YOU CALL YOUR SCRIPT AND COPY THE OUTPUT LOCALLY (so that it can get back into WDL-space)
                    submit_cmd = f'''gcloud dataproc jobs submit pyspark {script_path} \
                    --cluster={cluster_name} --project ~{gcs_project} --region=~{region} --account {account} --driver-log-levels root=WARN -- \
                    --input_filtered_vds_url ~{input_filtered_vds_url} \
                    --output_final_vcf_url ~{output_final_vcf_url} \
                    --output_final_vcf_header_url ~{output_final_vcf_header_url} \
                    --output_report_url ~{output_report_url} \
                    --contig ~{contig} \
                    --bed_url ~{bed_url} \
                    --temp_bucket gs://{cluster_temp_bucket}/{cluster_name}'''

                    print("Running: " + submit_cmd)
                    f = os.popen(submit_cmd)
                    f.read()
                    if (f.close() != None):
                        raise Exception("Failed to submit cluster job sucessfully")
                    ###########

                    break

        except Exception as e:
            print(e)
            raise
        finally:
            print(f'Stopping cluster: {cluster_name}')
            os.popen("gcloud dataproc clusters delete --project {} --region {} --account {} {}".format("~{gcs_project}", "~{region}", account, cluster_name)).read()

        EOF

        echo "Complete"
    >>>

    output {
        String final_vcf = output_final_vcf_url
        String final_vcf_header = output_final_vcf_header_url
        String report = output_report_url
    }

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}
