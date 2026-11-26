version 1.0
#
# AoU + AnVIL merge, filter and QC: from two VDSs and a BED file to one VCF per contig.
#
# Hail steps, each on its own Dataproc cluster:
#   Step 1  AoU site statistics   (every site; no filtering)   -> Hail Table + TSV     } in parallel
#   Step 2  AnVIL site statistics (every site; no filtering)   -> Hail Table + TSV     }
#   Step 3  apply the filters, restrict both VDSs to the kept sites, rename AnVIL samples, combine
#   Step 4  densify the combined VDS and export the final VCF + tabix index
# Optional:
#   Step 4 again on the AoU-only filtered VDS (export_aou_only): an AoU-only VCF at the same sites
#   EvalCallset + NormalizeEvalCallset (eval_sample_names): unfiltered genotypes of the samples
#     that have truth, for scoring filter sets offline with filter_sweep.py
#
# Typical use (see RUNBOOK.md):
#   A. run_site_stats = true,  run_panel = false, eval_sample_names set  -> statistics + evaluation callset
#   B. choose filters offline with filter_sweep.py
#   C. run_site_stats = false, run_panel = true, filter_spec / filter_tag set -> the panel, reusing A
#
# Filters (filter_spec): a comma-separated list of key=value with no spaces. "none" gives the
# defaults: AVSAD >= 12, call rate >= 0.90, mean GQ >= 30, AoU max alt AC >= 2, AoU LowQual /
# NO_HQ_GENOTYPES flags, AnVIL no-high-quality-genotype, and Hail's excess-het test at
# p < 3.4e-6 (per alt allele, on diploid calls; XX samples in chrX non-PAR, all elsewhere).
#   keys     <cohort>.<name> or <cohort>.<region>.<name>
#   cohort   aou | anvil | both          region   par | nonpar (omit for everywhere)
#   name     min_avsad | min_call_rate | min_gq | min_ac (aou) | excess_het_p | excess_het_f |
#            max_xy_het_frac | xy_het_min_carriers
#   also     aou.excess_het=hail|flag    (flag: AoU's precomputed ExcessHet flag outside non-PAR)
#   example  aou.par.min_call_rate=0.8,anvil.nonpar.min_avsad=6
# For autosomes, aou.excess_het=flag reproduces the original Step 1.
#
# Important notes:
#   - Every step returns its outputs as strings holding fixed bucket paths. A call-cache hit returns
#     the path without checking what is there. Steps 1 and 2 always write to the same names under
#     output_bucket_path, which is what lets a later run reuse them (run_site_stats = false).
#     Give every filter set its own filter_tag: Step 3's combiner refuses to run if a finished VDS
#     already exists at its output path.
#   - Cluster defaults are 2 workers + 8 preemptible workers of n2d-highmem-8. If you change
#     worker_machine_type, check executor_memory / executor_cores still fit it.
#   - All bed files referenced in this WDL are UCSC bed files (as opposed to PLINK bed files).
#   - This has not been tested on any reference other than hg38.
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
    input {

        ## ANALYSIS PARAMETERS
        #  ie, parameters that go to the Hail python code (submission scripts below)

        # Input files/paths
        String aou_vds_url
        String anvil_vds_url
        String bed_url

        # Sample sex for each cohort: GCS urls of two-column (sample, sex) tables, tab-, comma- or
        #  space-separated; XX / F / FEMALE and XY / M / MALE are recognised. The AnVIL table must use
        #  the sample IDs as they are in the AnVIL VDS (before the pipeline appends "_C").
        #  Only read when the contig is chrX. Used by Steps 1 and 2, and by Step 4 for masking.
        String aou_sex_tsv_url
        String anvil_sex_tsv_url

        # Contigs of interest.  If a contig is present in the bed file, but not in this list, the contig will be ignored.
        #   In other words, this is a contig level intersection with the bed file.
        #     This list of contigs that must be present in the reference.  Each contig will be processed separately (shard)
        # This list should be ordered.  Eg, ["chr21", "chr22"]
        Array[String] contigs

        # String used in construction of cluster names
        #  Cannot contain any special characters, ie, characters must be alphanumeric or "_"
        String prefix

        ## WHAT TO RUN
        # Steps 1 and 2. Set to false to reuse the statistics a previous run wrote under output_bucket_path.
        Boolean run_site_stats = true
        # Steps 3 and 4. Set to false to stop after the statistics (and the evaluation callset).
        Boolean run_panel = true
        # Filter thresholds for Step 3 ("none" = the defaults); see the notes at the top.
        String filter_spec = "none"
        # Short label for this filter set, inserted into the Step 3 and Step 4 output names so that
        #  several filter sets can share one output_bucket_path. Alphanumeric, "-" or "_". "" = no label,
        #  which keeps the final VCF name of the old pipeline.
        String filter_tag = ""
        # chrX only, Step 4: write heterozygous calls of XY samples in non-PAR as missing (./.), so that
        #  Beagle imputes them. A male has one X there, so such a call cannot be right; AnVIL males are
        #  diploid-encoded and carry many of them. Set to false for the genotypes exactly as called.
        Boolean mask_xy_nonpar_hets = true
        # Also export an AoU-only VCF at the same kept sites (for measuring what adding AnVIL buys)
        Boolean export_aou_only = false
        # AoU samples that have truth. If given, their unfiltered genotypes are written as a small
        #  normalised VCF for filter_sweep.py; needs ref_fasta and ref_fasta_index.
        Array[String] eval_sample_names = []
        File? ref_fasta
        File? ref_fasta_index

        ## CLUSTER PARAMETERS

        # The Google project ID information is necessary when spinning up dataproc.
        #  ie, terra-<hex>
        # This must match the workspace that this workflow is being run from.
        #     eg, "terra-491d5f31"
        String gcs_project

        # Set to 'subnetwork' if running in Terra Cromwell
        String gcs_subnetwork_name = 'subnetwork'

        # The scripts that are run on the clusters
        File submission_script_site_stats     # Step1_2_Site_Stats.py (Steps 1 and 2)
        File submission_script_step3          # Step3_Filter_Rename_Combine.py
        File submission_script_step4          # Step4_Densify_Export.py
        File submission_script_eval_callset   # Eval_Callset.py

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
    String out = sub(output_bucket_path, "/$", "") + "/"
    String tag_part = if filter_tag == "" then "" else filter_tag + "."

    scatter (contig in contigs) {

        # Fixed names, so that a later run can find the statistics without re-running Steps 1 and 2
        String aou_stats_ht    = out + contig + ".step1.aou.site_stats.ht"
        String aou_stats_tsv   = out + contig + ".step1.aou.site_stats.tsv.bgz"
        String anvil_stats_ht  = out + contig + ".step2.anvil.site_stats.ht"
        String anvil_stats_tsv = out + contig + ".step2.anvil.site_stats.tsv.bgz"

        if (run_site_stats) {
            call Step1 {
                input:
                    input_vds_url = aou_vds_url,
                    sex_tsv_url = aou_sex_tsv_url,
                    output_stats_ht_url = aou_stats_ht,
                    output_stats_tsv_url = aou_stats_tsv,
                    output_report_url = out + contig + ".step1.report",
                    contig = contig,
                    bed_url = bed_url,
                    prefix = prefix,
                    gcs_project = gcs_project,
                    gcs_subnetwork_name = gcs_subnetwork_name,
                    submission_script = submission_script_site_stats,
                    hail_docker = hail_docker,
                    region = region
            }

            call Step2 {
                input:
                    input_vds_url = anvil_vds_url,
                    sex_tsv_url = anvil_sex_tsv_url,
                    output_stats_ht_url = anvil_stats_ht,
                    output_stats_tsv_url = anvil_stats_tsv,
                    output_report_url = out + contig + ".step2.report",
                    contig = contig,
                    bed_url = bed_url,
                    prefix = prefix,
                    gcs_project = gcs_project,
                    gcs_subnetwork_name = gcs_subnetwork_name,
                    submission_script = submission_script_site_stats,
                    hail_docker = hail_docker,
                    region = region
            }
        }

        if (run_panel) {
            # Step 3 waits for Steps 1 and 2 when they run here, and otherwise reads what a previous run wrote
            call Step3 {
                input:
                    input_aou_vds_url = aou_vds_url,
                    input_anvil_vds_url = anvil_vds_url,
                    input_aou_stats_ht_url = select_first([Step1.stats_ht, aou_stats_ht]),
                    input_anvil_stats_ht_url = select_first([Step2.stats_ht, anvil_stats_ht]),
                    filter_spec = filter_spec,
                    tag_part = tag_part,
                    contig = contig,
                    bed_url = bed_url,
                    prefix = prefix,
                    gcs_project = gcs_project,
                    gcs_subnetwork_name = gcs_subnetwork_name,
                    submission_script = submission_script_step3,
                    output_bucket = out,
                    hail_docker = hail_docker,
                    region = region
            }

            call Step4 {
                input:
                    input_vds_url = Step3.combined_vds,
                    mask_xy_nonpar_hets = mask_xy_nonpar_hets,
                    aou_sex_tsv_url = aou_sex_tsv_url,
                    anvil_sex_tsv_url = anvil_sex_tsv_url,
                    tag_part = tag_part,
                    contig = contig,
                    prefix = prefix,
                    gcs_project = gcs_project,
                    gcs_subnetwork_name = gcs_subnetwork_name,
                    submission_script = submission_script_step4,
                    output_bucket = out,
                    hail_docker = hail_docker,
                    region = region
            }

            if (export_aou_only) {
                call Step4 as Step4AouOnly {
                    input:
                        input_vds_url = Step3.aou_filtered_vds,
                        mask_xy_nonpar_hets = mask_xy_nonpar_hets,
                        aou_sex_tsv_url = aou_sex_tsv_url,
                        anvil_sex_tsv_url = anvil_sex_tsv_url,
                        name_prefix = "aou-only",
                        tag_part = tag_part,
                        contig = contig,
                        prefix = prefix,
                        gcs_project = gcs_project,
                        gcs_subnetwork_name = gcs_subnetwork_name,
                        submission_script = submission_script_step4,
                        output_bucket = out,
                        hail_docker = hail_docker,
                        region = region
                }
            }
        }

        if (length(eval_sample_names) > 0) {
            call EvalCallset {
                input:
                    input_aou_vds_url = aou_vds_url,
                    sample_names = eval_sample_names,
                    contig = contig,
                    bed_url = bed_url,
                    prefix = prefix,
                    gcs_project = gcs_project,
                    gcs_subnetwork_name = gcs_subnetwork_name,
                    submission_script = submission_script_eval_callset,
                    output_bucket = out,
                    hail_docker = hail_docker,
                    region = region
            }

            call NormalizeEvalCallset {
                input:
                    vcf = EvalCallset.vcf,
                    ref_fasta = select_first([ref_fasta]),
                    ref_fasta_index = select_first([ref_fasta_index]),
                    output_basename = "eval_callset.unfiltered.norm." + contig
            }
        }
    }

    output {
        # Steps 1 and 2: where the site statistics are (written by this run, or reused)
        Array[String] aou_site_stats_hts = aou_stats_ht
        Array[String] aou_site_stats_tsvs = aou_stats_tsv
        Array[String] anvil_site_stats_hts = anvil_stats_ht
        Array[String] anvil_site_stats_tsvs = anvil_stats_tsv

        # Step 3
        Array[String?] aou_pass_hts = Step3.aou_pass_ht
        Array[String?] anvil_exclusion_hts = Step3.anvil_exclusion_ht
        Array[String?] keep_hts = Step3.keep_ht
        Array[String?] aou_filtered_vdss = Step3.aou_filtered_vds
        Array[String?] anvil_filtered_vdss = Step3.anvil_filtered_vds
        Array[String?] combined_vdss = Step3.combined_vds

        # Step 4 (was Step 5b)
        Array[String?] final_vcf = Step4.final_vcf
        Array[String?] final_vcf_index = Step4.final_vcf_index
        Array[String?] aou_only_vcf = Step4AouOnly.final_vcf
        Array[String?] aou_only_vcf_index = Step4AouOnly.final_vcf_index

        # Evaluation callset: unfiltered, normalised genotypes of eval_sample_names
        Array[File?] eval_callset_vcf = NormalizeEvalCallset.output_vcf
        Array[File?] eval_callset_vcf_tbi = NormalizeEvalCallset.output_vcf_tbi
    }
}

task Step1 {
    input {
        # You must treat a VDS as a String, since it is a directory and not a single file
        String input_vds_url
        String bed_url
        # Two-column (sample, sex) table for this cohort, sample IDs as in the VDS
        String sex_tsv_url

        # Where the outputs go (set by the workflow, so that a later run can reuse them)
        String output_stats_ht_url
        String output_stats_tsv_url
        String output_report_url
        # How far before each interval to look for reference blocks when a VDS has no stored
        #  maximum reference-block length (the AnVIL VDS does not).
        Int ref_pad_bp = 1000000

        # Cannot make this localization_optional
        File submission_script

        # contig must be in the reference
        String contig
        String prefix

        # dataproc params
        String gcs_project
        String region = "us-central1"
        String master_machine_type = "n1-highmem-32"
        Float master_memory_fraction = 0.8
        String worker_machine_type = "n2d-highmem-8"
        Int num_workers = 2
        Int num_preemptible_workers = 8
        Int time_to_live_minutes = 2880 # two days
        # Spark executor sizing, passed to the script. Must fit the worker machine type.
        String executor_memory = "22g"
        String executor_cores = "4"
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
                    --cohort aou \
                    --input_vds_url ~{input_vds_url} \
                    --sex_tsv_url ~{sex_tsv_url} \
                    --output_stats_ht_url ~{output_stats_ht_url} \
                    --output_stats_tsv_url ~{output_stats_tsv_url} \
                    --output_report_url ~{output_report_url} \
                    --bed_url ~{bed_url} \
                    --ref_pad_bp ~{ref_pad_bp} \
                    --contig ~{contig} \
                    --executor_memory ~{executor_memory} \
                    --executor_cores ~{executor_cores} \
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
        String stats_ht = output_stats_ht_url
        String stats_tsv = output_stats_tsv_url
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
        String input_vds_url
        String bed_url
        # Two-column (sample, sex) table for this cohort, sample IDs as in the VDS
        String sex_tsv_url

        # Where the outputs go (set by the workflow, so that a later run can reuse them)
        String output_stats_ht_url
        String output_stats_tsv_url
        String output_report_url
        # How far before each interval to look for reference blocks when a VDS has no stored
        #  maximum reference-block length (the AnVIL VDS does not).
        Int ref_pad_bp = 1000000

        # Cannot make this localization_optional
        File submission_script

        # contig must be in the reference
        String contig
        String prefix

        # dataproc params
        String gcs_project
        String region = "us-central1"
        String master_machine_type = "n1-highmem-32"
        Float master_memory_fraction = 0.8
        String worker_machine_type = "n2d-highmem-8"
        Int num_workers = 2
        Int num_preemptible_workers = 8
        Int time_to_live_minutes = 2880 # two days
        # Spark executor sizing, passed to the script. Must fit the worker machine type.
        String executor_memory = "22g"
        String executor_cores = "4"
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
                    --cohort anvil \
                    --input_vds_url ~{input_vds_url} \
                    --sex_tsv_url ~{sex_tsv_url} \
                    --output_stats_ht_url ~{output_stats_ht_url} \
                    --output_stats_tsv_url ~{output_stats_tsv_url} \
                    --output_report_url ~{output_report_url} \
                    --bed_url ~{bed_url} \
                    --ref_pad_bp ~{ref_pad_bp} \
                    --contig ~{contig} \
                    --executor_memory ~{executor_memory} \
                    --executor_cores ~{executor_cores} \
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
        String stats_ht = output_stats_ht_url
        String stats_tsv = output_stats_tsv_url
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

        String input_aou_stats_ht_url        # output from Step 1
        String input_anvil_stats_ht_url      # output from Step 2
        String bed_url

        # Filter thresholds ("none" = the defaults); see the notes at the top of this file
        String filter_spec
        # "" or "<filter_tag>.", inserted into the output names
        String tag_part
        String output_bucket
        # How far before each interval to look for reference blocks when a VDS has no stored
        #  maximum reference-block length (the AnVIL VDS does not).
        Int ref_pad_bp = 1000000

        # Cannot make this localization_optional
        File submission_script

        # contig must be in the reference
        String contig
        String prefix

        # dataproc params
        String gcs_project
        String region = "us-central1"
        String master_machine_type = "n1-highmem-32"
        Float master_memory_fraction = 0.8
        String worker_machine_type = "n2d-highmem-8"
        Int num_workers = 2
        Int num_preemptible_workers = 8
        Int time_to_live_minutes = 10080 # seven days
        # Spark executor sizing, passed to the script. Must fit the worker machine type.
        String executor_memory = "22g"
        String executor_cores = "4"
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
    String output_aou_pass_ht_url        = output_bucket + contig + ".step3." + tag_part + "aou.pass.ht"
    String output_anvil_exclusion_ht_url = output_bucket + contig + ".step3." + tag_part + "anvil.exclusion.ht"
    String output_keep_ht_url            = output_bucket + contig + ".step3." + tag_part + "keep.ht"
    String output_aou_filtered_vds_url   = output_bucket + contig + ".step3." + tag_part + "aou.filtered.vds"
    String output_anvil_filtered_vds_url = output_bucket + contig + ".step3." + tag_part + "anvil.filtered.renamed.vds"
    String output_combined_vds_url       = output_bucket + contig + ".step3." + tag_part + "aou-anvil.combined.vds"
    String output_report_url             = output_bucket + contig + ".step3." + tag_part + "report"

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
                    --filter_spec ~{filter_spec} \
                    --input_aou_stats_ht_url ~{input_aou_stats_ht_url} \
                    --input_anvil_stats_ht_url ~{input_anvil_stats_ht_url} \
                    --input_aou_vds_url ~{input_aou_vds_url} \
                    --input_anvil_vds_url ~{input_anvil_vds_url} \
                    --output_aou_pass_ht_url ~{output_aou_pass_ht_url} \
                    --output_anvil_exclusion_ht_url ~{output_anvil_exclusion_ht_url} \
                    --output_keep_ht_url ~{output_keep_ht_url} \
                    --output_aou_filtered_vds_url ~{output_aou_filtered_vds_url} \
                    --output_anvil_filtered_vds_url ~{output_anvil_filtered_vds_url} \
                    --output_combined_vds_url ~{output_combined_vds_url} \
                    --output_report_url ~{output_report_url} \
                    --bed_url ~{bed_url} \
                    --ref_pad_bp ~{ref_pad_bp} \
                    --contig ~{contig} \
                    --executor_memory ~{executor_memory} \
                    --executor_cores ~{executor_cores} \
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
        String aou_pass_ht = output_aou_pass_ht_url
        String anvil_exclusion_ht = output_anvil_exclusion_ht_url
        String keep_ht = output_keep_ht_url
        String aou_filtered_vds = output_aou_filtered_vds_url
        String anvil_filtered_vds = output_anvil_filtered_vds_url
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

task Step4 {
    input {
        # You must treat a VDS as a String, since it is a directory and not a single file
        String input_vds_url               # combined VDS (or the filtered AoU-only VDS) from Step 3

        # chrX only: write heterozygous calls of XY samples in non-PAR as missing, so Beagle imputes them
        Boolean mask_xy_nonpar_hets
        # Two-column (sample, sex) tables, as for Steps 1 and 2 (AnVIL IDs as in the AnVIL VDS)
        String aou_sex_tsv_url
        String anvil_sex_tsv_url

        # First part of the output file names: "aou-anvil" for the panel, "aou-only" for the comparison arm
        String name_prefix = "aou-anvil"
        # "" or "<filter_tag>.", inserted into the output names
        String tag_part
        String output_bucket

        # Cannot make this localization_optional
        File submission_script

        # contig must be in the reference
        String contig
        String prefix

        # dataproc params
        String gcs_project
        String region = "us-central1"
        String master_machine_type = "n1-highmem-64"
        Float master_memory_fraction = 0.8
        String worker_machine_type = "n2d-highmem-8"
        Int num_workers = 2
        Int num_preemptible_workers = 8
        Int time_to_live_minutes = 2880 # two days
        # Spark executor sizing, passed to the script. Must fit the worker machine type.
        String executor_memory = "45g"
        String executor_cores = "8"
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
    # With the default name_prefix and no filter_tag the file names are those of the old Step 5b
    String output_final_vcf_url        = output_bucket + name_prefix + ".final.step5b." + tag_part + contig + ".vcf.bgz" # imputation wdl expects contig right before prefix
    String output_final_vcf_header_url = output_bucket + name_prefix + ".final." + tag_part + contig + ".step5b.vcf.header"
    String output_report_url           = output_bucket + contig + ".step4." + tag_part + name_prefix + ".report"

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
                    --input_combined_vds_url ~{input_vds_url} \
                    --output_final_vcf_url ~{output_final_vcf_url} \
                    --output_final_vcf_header_url ~{output_final_vcf_header_url} \
                    --output_report_url ~{output_report_url} \
                    --mask_xy_nonpar_hets ~{mask_xy_nonpar_hets} \
                    --aou_sex_tsv_url ~{aou_sex_tsv_url} \
                    --anvil_sex_tsv_url ~{anvil_sex_tsv_url} \
                    --contig ~{contig} \
                    --executor_memory ~{executor_memory} \
                    --executor_cores ~{executor_cores} \
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
        String final_vcf_index = output_final_vcf_url + ".tbi"
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

task EvalCallset {
    input {
        # You must treat a VDS as a String, since it is a directory and not a single file
        String input_aou_vds_url
        String bed_url
        Array[String] sample_names
        String output_bucket
        # How far before each interval to look for reference blocks when a VDS has no stored
        #  maximum reference-block length (the AnVIL VDS does not).
        Int ref_pad_bp = 1000000

        # Cannot make this localization_optional
        File submission_script

        # contig must be in the reference
        String contig
        String prefix

        # dataproc params
        String gcs_project
        String region = "us-central1"
        String master_machine_type = "n1-highmem-32"
        Float master_memory_fraction = 0.8
        String worker_machine_type = "n2d-highmem-8"
        Int num_workers = 2
        Int num_preemptible_workers = 8
        Int time_to_live_minutes = 2880 # two days
        # Spark executor sizing, passed to the script. Must fit the worker machine type.
        String executor_memory = "22g"
        String executor_cores = "4"
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
    String output_vcf_url    = output_bucket + contig + ".eval_callset.unfiltered.vcf.bgz"
    String output_report_url = output_bucket + contig + ".eval_callset.report"

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
        cluster_name = f'~{prefix}-~{contig}-hail-evalset-{str(uuid.uuid4())[0:13]}'.lower()

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
                    --samples ~{sep="," sample_names} \
                    --output_vcf_url ~{output_vcf_url} \
                    --output_report_url ~{output_report_url} \
                    --bed_url ~{bed_url} \
                    --ref_pad_bp ~{ref_pad_bp} \
                    --contig ~{contig} \
                    --executor_memory ~{executor_memory} \
                    --executor_cores ~{executor_cores} \
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
        String vcf = output_vcf_url
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

task NormalizeEvalCallset {
    input {
        File vcf
        File ref_fasta
        File ref_fasta_index
        String output_basename

        Int cpu = 4
        Int memory_gb = 16
        Int disk_size_gb = ceil(size(vcf, "GB") * 6 + size(ref_fasta, "GB")) + 20
        Int preemptible_tries = 2
    }

    command <<<
        set -euo pipefail

        # bcftools needs the .fai next to the FASTA
        ln -s ~{ref_fasta} ref.fa
        ln -s ~{ref_fasta_index} ref.fa.fai

        # Same normalisation as BeaglePhasing.wdl applies to the panel:
        #   1. split multiallelic sites
        #   2. left-align against the reference, then drop exact duplicates
        #   3. drop records that no evaluation sample carries, compress and index
        bcftools norm -m - -Ou ~{vcf} \
            | bcftools norm -f ref.fa --rm-dup exact -Ou \
            | bcftools view --min-ac 1 -Oz --threads ~{cpu} --write-index=tbi -o ~{output_basename}.vcf.gz
    >>>

    output {
        File output_vcf = "~{output_basename}.vcf.gz"
        File output_vcf_tbi = "~{output_basename}.vcf.gz.tbi"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: preemptible_tries
    }
}
