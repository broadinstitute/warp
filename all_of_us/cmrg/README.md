# Challenging Medically Relevant Genes Variant Calling Workflows (WDL)

This directory contains WDL workflows used for **Challenging Medically Relevant Genes (CMRG) analysis** variant calling and preparation steps. In broad strokes, these workflows support:

1. Subsetting reads to CMRG-relevant regions
2. Remapping reads in problematic (“false duplication”) regions using **FixItFelix** against a masked GRCh38 reference
3. Calling variants with **GATK HaplotypeCaller** on corrected alignments
4. (Optional downstream) reblocking GVCFs for efficient storage / joint-calling ingestion
5. (External) joint variant calling and extraction via the **GATK Variant Store (GVS)** pipeline

---

## Workflow 1: `FixItFelixAndVariantCall`

### Description

End-to-end workflow that:

1. **Subsets** an input CRAM to intervals spanning “true locations” and “false duplications”
2. Runs **FixItFelix** to extract reads in these regions, convert to FASTQ, and **remap** to a **masked GRCh38 reference**
3. Runs **GATK HaplotypeCaller** on the corrected BAM in the *true* location intervals and outputs a VCF (or optional gVCF)

### Required Inputs

* `cram_file` (File): Input CRAM (or BAM; workflow logic supports either)
* `cram_file_index` (File): CRAI (or BAI)
* `original_ref_fasta` / `original_ref_fasta_index` / `original_ref_dict` (File): Reference used for the original CRAM/BAM alignment

Defaults are provided in the WDL for:

* Masked GRCh38 reference + indexes (`masked_ref_*`)
* Interval BEDs:

  * `true_locations_intervals`
  * `false_duplications_intervals`
  * `combined_true_false_intervals`

Optional:

* `generate_gvcf` (Boolean, default `false`): if `true`, emits gVCF output

### Tasks & Software

#### Task: `subset_cram`

* **Tool:** GATK `PrintReads`
* **Docker:** `us.gcr.io/broad-gatk/gatk:4.4.0.0`
* **Purpose:** Extract only reads overlapping the combined true/false interval BED and output a BAM

#### Task: `FixItFelix`

* **Tools:** `samtools`, `bwa mem`, helper utilities within container
* **Docker:** `gcr.io/broad-dsde-methods/fixitfelix:1`
* **Purpose:** For the provided intervals:

  * extract properly paired reads
  * convert to FASTQ
  * remap to the **masked** GRCh38 reference
  * output a remapped, coordinate-sorted, indexed BAM

> Note: the task has an `old_ref` input (default empty). The script enforces that CRAM inputs require an old reference. In this repo’s integrated workflow, CRAMs are first converted to BAM via `subset_cram`, so `old_ref` is typically not needed.

#### Task: `call_variants`

* **Tool:** GATK `HaplotypeCaller` (with `--dragen-mode`)
* **Docker:** `us.gcr.io/broad-gatk/gatk:4.4.0.0`
* **Purpose:**

  * Call variants on the remapped BAM in the **true location** intervals
  * If producing a VCF (not gVCF), run `SelectVariants --exclude-non-variants` to drop non-variant sites

### Outputs

* `output_vcf` (File): `*.vcf.gz` (filtered) or `*.g.vcf.gz` (if `generate_gvcf=true`)
* `output_vcf_index` (File): `*.tbi`
* `output_pipeline_version` (String): set in workflow (`aou_9.0.1`)

---

## Workflow 2: `GvsAoUReblockGvcf`

### Description

Reblocks an existing gVCF (commonly sourced from DRC buckets) using **GATK ReblockGVCF**, and optionally copies the reblocked output to an AoU research bucket location based on `site_id`.

This workflow intentionally combines localization + delocalization into a single task for cost efficiency.

### Required Inputs

* `gvcf` (String): GCS path to `*.g.vcf.gz`
* `gvcf_index` (String?, optional): defaults to `gvcf + ".tbi"`
* `ref_fasta` / `ref_fasta_index` / `ref_dict` (File): Reference inputs for GATK
* `docker_image` (String, default `us.gcr.io/broad-gatk/gatk:4.2.6.1`)

Optional routing / billing:

* `site_id` (String?, optional): one of `bi`, `bcm`, `uw` (case-normalized in WDL)
* `requester_pays_project` (String?, optional): passed to GATK as requester pays project

### Tasks & Software

#### Task: `ReblockAndCopy`

* **Tool:** GATK `ReblockGVCF`
* **Docker:** `docker_image` (default: `us.gcr.io/broad-gatk/gatk:4.2.6.1`)
* **Optional copy:** `gsutil -m cp` to the derived site bucket/prefix when `site_id` is provided

### Outputs

* `reblocked_gvcf` (File): destination path + `*.reblocked.g.vcf.gz`
* `reblocked_gvcf_index` (File): destination path + `*.tbi`

> If `site_id` is not provided, the intent (per comments) is to “drop locally” rather than copy; however, note the output section currently constructs `destination + output_gvcf_filename`, so callers should ensure outputs match the desired behavior for their execution environment.

---

## Downstream Joint Calling (External): `GvsJointVariantCalling` (GVS Store)

### Description (high-level only)

Downstream of per-sample gVCF creation/reblocking, joint calling is performed using the **GATK Variant Store (GVS)** workflow **outside this repository**. At a high level, it:

* Bulk ingests per-sample VCF/gVCF inputs into GVS tables
* Populates alternate allele information
* Creates a filter set (e.g., VETS/VQSR depending on configuration)
* Prepares and extracts cohort callsets (optionally merged/bgzip’d), producing final VCF outputs and manifests

The canonical workflow is maintained in the GATK repository. ([GitHub][1])

For convenience, the source location is:

```text
https://github.com/broadinstitute/gatk/blob/gvs_0.6.4/scripts/variantstore/wdl/GvsJointVariantCalling.wdl
```

---

## Suggested Usage Pattern

A common sequence is:

1. Run **FixItFelixAndVariantCall** to produce corrected-region VCFs (and/or gVCFs if enabled)
2. If needed for joint calling/storage, run **GvsAoUReblockGvcf** on gVCFs
3. Use the external **GvsJointVariantCalling** workflow to ingest and generate joint callsets

