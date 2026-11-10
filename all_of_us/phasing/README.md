# Phasing Pipelines
The following pipelines are used to phase WGS data.

## filter_and_qc_variants
#### Background

This WDL workflow processes genomic data by filtering and performing quality control (QC) on variants. 
It uses Hail and Google Cloud Dataproc to handle large-scale genomic datasets efficiently. The workflow is tailored 
for the hg38 reference genome. Each contig handled separately and requires a separate run of the workflow 

Key characteristics:
- Each chromosome is processed independently.
- Output VCFs are not merged.
- Designed for scalability (e.g., 500k AoU WGS VDS).
- Not compatible with requester-pays buckets in Cromwell.
- Input Data: Accepts a Variant Dataset (VDS) and a UCSC BED file to define genomic regions of interest.

#### Inputs
Analysis Parameters:
- `String input_aou_vds_url` – URL of the input VDS
- `File submission_script` – Python script executed on the Dataproc cluster; defines filter and QC steps
- `String output_bucket` – GCS bucket for storing outputs
- `String contig` – Contig to process (e.g., "chr21")
- `String prefix` – Prefix for output filenames

Cluster Parameters:
- `String gcs_project` – Google Cloud project ID
- `String region` – Google Cloud region (default: "us-central1")
- `String master_machine_type` – Machine type for the master node (default: "n1-highmem-32")
- `Float master_memory_fraction` – Memory fraction for the master node (default: 0.8)
- `String worker_machine_type` – Machine type for worker nodes (default: "n1-highmem-4")
- `Int num_workers` – Number of worker nodes (default: 2)
- `Int num_preemptible_workers` – Number of preemptible worker nodes (default: 50)
- `Int time_to_live_minutes` – Time to live for the cluster (default: 2880 minutes)
- `String gcs_subnetwork_name` – Subnetwork name for Dataproc

Runtime Parameters:
- `RuntimeAttr? runtime_attr_override` – Optional runtime attributes for the job
- `String hail_docker` – Docker image with Hail and Google Cloud SDK

#### Step 1. FilterAndQCVariants
- Input Validation: Ensures required inputs are provided and properly formatted.
- Cluster Setup: Spins up a Dataproc cluster with user-defined configurations.
- Data Processing:
  - Filters variants based on the provided BED file.
  - Performs QC on the filtered variants.
- Output Generation: Saves filtered VCFs, headers, and QC reports to the specified GCS bucket.
- Cluster Teardown: Deletes the Dataproc cluster after processing.

#### Outputs

- `String aou_vcf` – URL for single chromosome output VCF
- `String aou_vcf_header` – URL for single chromosome output VCF header
- `String report` – URL for the report file

## Beagle5Phasing
#### Background

This WDL workflow phases genomic variants using Beagle5. It is designed for use in the All of Us Local Ancestry
pipeline and operates on a single VCF at a time. A user-supplied shell script drives the Beagle5 execution,
allowing flexibility in phasing parameters.

Key characteristics:
- Each VCF is processed independently; outputs are not merged across chromosomes.
- Phasing is performed with a sliding-window approach controlled by the `window_markers` parameter.
- Designed for large-scale WGS datasets (e.g., AoU cohort).
- Input Data: Accepts a GCS-hosted VCF and a genetic map file.

#### Inputs
Analysis Parameters:
- `String vcf_gs` – GCS URL of the input VCF to phase
- `String map_gs` – GCS URL of the genetic map file
- `String out_prefix` – Prefix for output filenames
- `String out_dir_gs` – GCS directory for storing outputs
- `Int window_markers` – Window size in markers for phasing (default: 3500000)
- `String java_xmx` – Java max heap size for Beagle5 (default: "80g")
- `File run_beagle_sh` – Shell script executed to invoke Beagle5; defines phasing steps

Runtime Parameters:
- `Int runtime_cpu` – Number of CPUs (default: 32)
- `Int mem_gb` – Memory in GB (default: 120)
- `Int runtime_disk_gb` – Disk size in GB (default: 2000)
- `String runtime_disk_type` – Disk type, "HDD" or "SSD" (default: "HDD")

#### Step 1. PhaseWithBeagle
- Script Execution: Runs the user-supplied shell script with the provided phasing parameters.
- Phasing: Invokes Beagle5 to phase the input VCF using the supplied genetic map and window size.
- Output Generation: Saves the phased VCF and log file to the specified GCS output directory.

#### Outputs

- `String phased_vcf_gs` – GCS URL for the phased output VCF
- `String phased_log_gs` – GCS URL for the Beagle5 phasing log file
