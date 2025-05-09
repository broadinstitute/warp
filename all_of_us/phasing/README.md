# Phasing Pipelines
The following pipelines are used to phase WGS data.

## filter_and_qc_variants
#### Background

This WDL workflow processes genomic data by filtering and performing quality control (QC) on variants. 
It uses Hail and Google Cloud Dataproc to handle large-scale genomic datasets efficiently. The workflow is tailored 
for the hg38 reference genome and processes data in a sharded manner, with each contig handled separately.  

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



