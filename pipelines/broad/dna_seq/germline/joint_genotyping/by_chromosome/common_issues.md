# Common issues with JointGenotypingByChromosome

## Permissions
Since reblocking is done in FireCloud/Terra, the correct permissions must be added to the workspace in order for the workflow to read the input GVCFs
We are using the `broad-exomes-prodN`, where `N` is a number 1 through 10 inclusive, and the service accounts for each project must be added to the workspace as reader accounts.
The Cromwell execution account should also be added. The service accounts are listed below:
```text
420463522347-compute@developer.gserviceaccount.com
553968312475-compute@developer.gserviceaccount.com
22318392461-compute@developer.gserviceaccount.com
5753440872-compute@developer.gserviceaccount.com
1085246091060-compute@developer.gserviceaccount.com
572127813794-compute@developer.gserviceaccount.com
124548794588-compute@developer.gserviceaccount.com
561786753150-compute@developer.gserviceaccount.com
1059520358416-compute@developer.gserviceaccount.com
1052986742836-compute@developer.gserviceaccount.com
477026979441-compute@developer.gserviceaccount.com
```

## ImportGVCFs and GnarlyGenotyper failures
The most common failure for `ImportGVCFs` and `GnarlyGenotyper` is for the tools to run out of memory. The memory required for these tasks is the jvm memory plus memory needed for GenomicsDB to operate.
We have not seen the JVM run out of memory. When GenomicsDB runs out of memory, it does so silently. If shards for either of these tasks fails without any logging as to detailing why, upping the amount of memory 
in the `runtime {...}` block should clear up the problem. Make sure that the `command {...}` block isn't modified, or else call-caching will break for all other shards.

A rarer error is for GenomicsDB to just give up no matter how much memory is thrown at it. We've seen success with manually running the tasks and providing the data through the 
`workaround` mechanisms. The stderr of the task should give you the command to run. You'll need to have local versions of the reference fasta, reference fasta index, and reference dictionary handy to pass to the tools.
The ImportGVCFs workarounds are keyed on the `callset_name.splitted_interval_list_name`, while the GnarlyGenotyper workarounds are keyed on the `callset_name.top_level_shard_index.gnarly_shard_index`.
For GnarlyGenotyper, the annotationDBs and output VCFs need to be provided.
Once the data is generated, it can be provided to the WDL. Below is how we did this for the UKBB effort:
Once the data is generated, it can be provided to the WDL. Below is how we did this for the UKBB effort:
```json
{
  "JointGenotypingByChromosomePartOne.gnarly_workaround": {
  "ukbb_exome_joint_genotyping_batch_1_chr1.0.8": "gs://broad-pharma5-ukbb-inputs-broad/batch1/workarounds/ukbb_exome_joint_genotyping_batch_1_chr1.0.8.vcf.gz",
  "ukbb_exome_joint_genotyping_batch_1_chr2.38.3": "gs://broad-pharma5-ukbb-inputs-broad/batch1/workarounds/ukbb_exome_joint_genotyping_batch_1_chr2.38.3.vcf.gz",
  "ukbb_exome_joint_genotyping_batch_1_chr3.10.7": "gs://broad-pharma5-ukbb-inputs-broad/batch1/workarounds/ukbb_exome_joint_genotyping_batch_1_chr3.10.7.vcf.gz",
  "ukbb_exome_joint_genotyping_batch_1_chr8.46.4": "gs://broad-pharma5-ukbb-inputs-broad/batch1/workarounds/ukbb_exome_joint_genotyping_batch_1_chr8.46.4.vcf.gz",
  "ukbb_exome_joint_genotyping_batch_1_chr11.38.3": "gs://broad-pharma5-ukbb-inputs-broad/batch1/workarounds/ukbb_exome_joint_genotyping_batch_1_chr11.38.3.vcf.gz"
},
  "JointGenotypingByChromosomePartOne.gnarly_workaround_annotations": {
    "ukbb_exome_joint_genotyping_batch_1_chr1.0.8": "gs://broad-pharma5-ukbb-inputs-broad/batch1/workarounds/ukbb_exome_joint_genotyping_batch_1_chr1.0.8.annotations.vcf.gz",
    "ukbb_exome_joint_genotyping_batch_1_chr2.38.3": "gs://broad-pharma5-ukbb-inputs-broad/batch1/workarounds/ukbb_exome_joint_genotyping_batch_1_chr2.38.3.annotations.vcf.gz",
    "ukbb_exome_joint_genotyping_batch_1_chr3.10.7": "gs://broad-pharma5-ukbb-inputs-broad/batch1/workarounds/ukbb_exome_joint_genotyping_batch_1_chr3.10.7.annotations.vcf.gz",
    "ukbb_exome_joint_genotyping_batch_1_chr8.46.4": "gs://broad-pharma5-ukbb-inputs-broad/batch1/workarounds/ukbb_exome_joint_genotyping_batch_1_chr8.46.4.annotations.vcf.gz",
    "ukbb_exome_joint_genotyping_batch_1_chr11.38.3": "gs://broad-pharma5-ukbb-inputs-broad/batch1/workarounds/ukbb_exome_joint_genotyping_batch_1_chr11.38.3.annotations.vcf.gz"
  },
  "JointGenotypingByChromosomePartOne.gnarly_workaround_keys": "ukbb_exome_joint_genotyping_batch_1_chr1.0.8|ukbb_exome_joint_genotyping_batch_1_chr2.38.3|ukbb_exome_joint_genotyping_batch_1_chr3.10.7|ukbb_exome_joint_genotyping_batch_1_chr8.46.4|ukbb_exome_joint_genotyping_batch_1_chr11.38.3",
  "JointGenotypingByChromosomePartOne.import_workaround": {
    "ukbb_exome_joint_genotyping_batch_1_chr2.0001-scattered.interval_list": "gs://broad-pharma5-ukbb-inputs-broad/batch1/workarounds/chr2/genomicsdb.tar"
  },
  "JointGenotypingByChromosomePartOne.import_workaround_keys": "ukbb_exome_joint_genotyping_batch_1_chr2.0001-scattered.interval_list"
 }
```

## Diagnosing Errors
The `by_chromosome_client.sh` script can give a nice output summarizing the status of each task in the WDL. It also gives the fully qualified name for tasks, which is useful for getting errors.
Given the `monitor` output of:
```text
chr1    3477d560-c2aa-43ae-8b95-c7297a396ed7    Running broad-exomes-prod5
        JointGenotypingByChromosomePartOne.CheckSamplesUnique      0 Running, 1 Done, 0 Preempted 0 Failed
        JointGenotypingByChromosomePartOne.GetFingerprintingIntervalIndices        0 Running, 1 Done, 0 Preempted 0 Failed
        JointGenotypingByChromosomePartOne.GnarlyIntervalScatterDude       0 Running, 78 Done, 0 Preempted 0 Failed
        JointGenotypingByChromosomePartOne.ImportGVCFs     77 Running, 0 Done, 0 Preempted 1 Failed
        JointGenotypingByChromosomePartOne.SplitIntervalList       0 Running, 1 Done, 0 Preempted 0 Failed
```
If a ImportGVCFs shard failed, a useful `jq` way of getting the stderr is `./by_chromosome_client.sh metadata -c chr1 | jq -r '.calls["ExomeJointGenotypingByChromosomePartOne.ImportGVCFs"][] | select(.executionStatus=="Failed") | .stderr' | less`

Sub-scatters can be a little trickier. Given the `monitor -x` output:
```text
chr1	cff51705-eea2-4577-a2ce-5f2e6fee4aa5	Succeeded	broad-pharma5-compute1
	eJointGenotypingByChromosomePartOne.CheckSamplesUnique	0 Running, 1 Done, 0 Preempted 0 Failed
	eJointGenotypingByChromosomePartOne.GatherAnnotationDBVcf	0 Running, 1 Done, 0 Preempted 0 Failed
	eJointGenotypingByChromosomePartOne.GatherFingerprintingVcfs	0 Running, 1 Done, 0 Preempted 0 Failed
	eJointGenotypingByChromosomePartOne.GetFingerprintingIntervalIndices	0 Running, 1 Done, 0 Preempted 0 Failed
	eJointGenotypingByChromosomePartOne.GnarlyIntervalScatterDude	0 Running, 82 Done, 0 Preempted 0 Failed
	eJointGenotypingByChromosomePartOne.HardFilterAndMakeSitesOnlyVcf	0 Running, 82 Done, 0 Preempted 0 Failed
	eJointGenotypingByChromosomePartOne.ImportGVCFs	0 Running, 82 Done, 0 Preempted 0 Failed
	eJointGenotypingByChromosomePartOne.SitesOnlyGatherVcf	0 Running, 1 Done, 0 Preempted 0 Failed
	eJointGenotypingByChromosomePartOne.SplitIntervalList	0 Running, 1 Done, 0 Preempted 0 Failed
	eJointGenotypingByChromosomePartOne.TotallyAwesomeGatherAnnotationDBs	0 Running, 82 Done, 0 Preempted 0 Failed
	eJointGenotypingByChromosomePartOne.TotallyRadicalGatherVcfs	0 Running, 82 Done, 0 Preempted 0 Failed
	SubWorkflow ScatterAt137_16
		JointGenotypingByChromosomePartOne.GnarlyGenotyper	0 Running, 9 Done, 0 Preempted 1 Failed
		JointGenotypingByChromosomePartOne.GnarlyGenotyper	0 Running, 10 Done, 0 Preempted 0 Failed
		JointGenotypingByChromosomePartOne.GnarlyGenotyper	0 Running, 10 Done, 0 Preempted 0 Failed
```
Getting the stderr for a failed GnarlyGenotyper shard would look like `./by_chromosome_client.sh metadata -x -c chr1 | jq -r '.calls["ScatterAt137_16"][].subWorkflowMetadata.calls["ExomeJointGenotypingByChromosomePartOne.GnarlyGenotyper"][] | select(.executionStatus=="Failed") | .stderr'`


## CrossCheckFingerprints INCONCLUSIVE Samples
Sometimes samples make it into the callset that have low coverage, and these samples will show up as `INCONCLUSIVE` in the `crosscheck_metrics`. The `CrossCheckFingerprints` task will fail,
but that doesn't mean that the workflow is doomed. Low-coverage samples have a statistically tiny chance of contributing anything to the final results, so the most common fix is to allow the sample to be `INCONCLUSIVE`.
This should only be done after confirming with a scientific owner or advisor. 

To allow an `INCONCLUSIVE` sample through, add a `"ExomeJointGenotypingByChromosomePartTwo.CrossCheckFingerprintsScattered.expected_inconclusive_samples": []` parameter to the `JointGenotypingByChromosomePartTwo.input.json`.
Any sample alias put into the array will be allowed through as `INCONCLUSIVE`.


## Manually Generating "Workaround" Data
Generating the data is a manual process, most easily done in GCP. Create a VM in one of the execution projects used for the real run. It should be a pretty big machine if being used for multiple workaround data generation jobs. An 8-core machine usually works well. Make sure there's more than enough disk space; around 100 GiB should do.
When provisioning a VM, don't deploy a container image to it. Just use the default Debian image and bump up the specs. This will give you `gsutil`. Deploying a container image will not give you `gsutil`, which is a pain to install. It also does not behave like the Google documentation says it should. Using the default Debian image provides a more sane environment to work in.
Once provisioned, [use this guide to install docker](https://docs.docker.com/install/linux/docker-ce/debian/). Make sure to configure permissions in the [post-install instructions](https://docs.docker.com/install/linux/linux-postinstall/) too.

References will need to be localized to the machine. The reference fasta, reference fasta index, and reference dictionary can be copied from `gs://gcp-public-data--broad-references/hg38/v0/`

An example of how to generate data for a failed GnarlyGenotyper shard is below. All the information needed for variables can be found in the shard's `script` 

```bash

mkdir -p ~/${CHROMOSOME}/${SHARD_NUMBER} && cd ~/${CHROMOSOME}/${SHARD_NUMBER}

mkdir ~/references

gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai \
  gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict \
  ~/references/

gsutil cp ${CLOUD_GENOMICS_DB} ./
tar -xf genomicsdb.tar && rm -f genomicsdb.tar

cat <<EOF > genomicsdb/query.json
  {
    "scan_full": true,
    "workspace": "genomicsdb",
    "array": "genomicsdb_array",
    "vid_mapping_file": "genomicsdb/vidmap.json",
    "callset_mapping_file": "genomicsdb/callset.json",
    "reference_genome": "/references/Homo_sapiens_assembly38.fasta",
    "max_diploid_alt_alleles_that_can_be_genotyped": 6,
    "produce_GT_field": true
  }
EOF

docker run -d \
  -v $(pwd):/workingDir \
  -v ~/references/:/references \
  ${GNARLY_DOCKER_IMAGE} \
  gatk --java-options "-Xms5000m -Xmx6500m"\
  GnarlyGenotyper \
  -R /references/Homo_sapiens_assembly38.fasta \
  -O /workingDir/${OUTPUT_VCF_BASE_NAME}.vcf.gz \
  --output-database-name /workingDir/${OUTPUT_VCF_BASE_NAME}.annotations.vcf.gz \
  -D gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf \
  --only-output-calls-starting-in-intervals \
  -V gendb:///workingDir/genomicsdb \
  -L ${CLOUD_INTERVALS_PATH} \
  -stand-call-conf 10 \
  --merge-input-intervals
```

## Becoming `picard-prod` for Submitting Part Two and Copying Cloud Outputs
To copy things around in `broad-gotc-prod-storage` or `broad-exomes-prod-storage` after manually generating data, the `picard-prod` service account needs to be used.

First, copy the `picard-prod` pem to the cloud VM. This can be done by using the gcloud SSH command provided by google, but replacing `ssh` with `scp` and providing paths like normal `scp`.
The file to copy can be generated by running `vault read -format=json -field=data secret/dsde/gotc/prod/picard/picard-account.pem > picard-account.json`.
The command will look something like `gcloud compute --project "broad-exomes-prod1" scp --zone "{ZONE}" picard-account.json "{VM_NAME}":~/picard-account.json`

To authenticate as `picard-prod`, run `gcloud auth activate-service-account --key-file {COPIED_PICARD_ACCOUNT_JSON}`.
After this, you will be authenticated as `picard-prod` and and can copy manually generated data to the workaround location.

To authenticate as your personal user again, run `gcloud config set account ${BROAD_EMAIL_ADDRESS}`. 
To switch back to `picard-prod`, run `gcloud config set account picard-prod@broad-gotc-prod.iam.gserviceaccount.com`
