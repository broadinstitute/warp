### Joint Calling on hg38 with the WARP JointGenotyping Pipeline

This workspace walks you through **joint calling** a cohort of samples on the hg38 (GRCh38) human reference, using the Broad Data Sciences Platform's cloud-optimized [JointGenotyping pipeline](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/dna_seq/germline/joint_genotyping/JointGenotyping.wdl). It is written for users who are **new to Terra and to joint calling**, and it takes you from per-sample GVCF files all the way to a filtered, multi-sample VCF that is ready to hand off to [**seqr**](https://seqr.broadinstitute.org/) for rare-disease variant analysis.

By the end of this tutorial you will be able to:

1. (Optional) **Reblock** your per-sample GVCF files to make joint calling faster and cheaper.
2. Build a **sample map** that tells the pipeline where your GVCFs live.
3. Run the **JointGenotyping** workflow in Terra to produce a filtered, joint-called VCF.
4. Understand how to **hand the joint VCF off to seqr**.

> **New to the terminology?** A **GVCF** ("genomic VCF") is a per-sample file that records evidence at *every* position in the genome, not just the variant sites. **Joint calling** combines many GVCFs so that variants are genotyped *together* across all samples — this improves sensitivity and gives you a single multi-sample VCF where every sample has a genotype at every site.

---

## Before You Begin

**What you need:**

- A **Terra account** and access to this workspace. If you are new to Terra, start with the [Terra Getting Started guide](https://support.terra.bio/hc/en-us/sections/23504885621787-Getting-Started-GCP).
- Your samples as **per-sample GVCF files** (plus their `.tbi` index files), called with GATK HaplotypeCaller in GVCF mode and aligned to **hg38 / GRCh38**. If you generated your samples with the [Whole Genome Germline Single Sample](https://broadinstitute.github.io/warp/docs/Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/README) pipeline, its GVCF output is exactly what you need (and it is already reblocked — see below).
- A **minimum of 50 samples.** The pipeline is designed for cohorts and will stop if fewer than ~50 unique samples are supplied. **Gene panels are not supported.**

**Reference data:** JointGenotyping needs a number of hg38 reference and resource files (reference FASTA, dbSNP, HapMap, Omni, 1000 Genomes, Mills, Axiom, interval lists, and a haplotype database for fingerprinting). **These are already configured for you in this workspace** — you do not need to supply them yourself.

---

## Workflow Overview

```
per-sample GVCFs
      │
      ▼
(0) ReblockGVCF        ← optional but HIGHLY recommended (faster + cheaper joint calling)
      │  reblocked GVCFs (*.rb.g.vcf.gz)
      ▼
(1) Generate a sample map   ← a 2-column list: sample name → GVCF path
      │  sample_map.tsv
      ▼
(2) JointGenotyping    ← joint calling + VQSR filtering (+ QC metrics, fingerprinting)
      │  <callset_name>.vcf.gz  (filtered, multi-sample)
      ▼
(3) Hand off to seqr
```

---

## Step 0 — (Recommended) Reblock your GVCFs

**Why reblock?** Reblocking condenses the reference (non-variant) blocks in a GVCF, producing a much smaller file **without losing variant information**. Smaller GVCFs mean joint calling runs **faster and costs less**, especially as your cohort grows. We **highly recommend** reblocking before joint calling.

> **You may be able to skip this step.** If your GVCFs already end in `.rb.g.vcf.gz`, or they came from the WARP Whole Genome Germline Single Sample pipeline (which reblocks by default), they are already reblocked — go straight to Step 1(below).

The [ReblockGVCF workflow](https://broadinstitute.github.io/warp/docs/Pipelines/ReblockGVCF_Pipeline/README) runs on **one sample at a time**, so in Terra you run it over a data table of samples and it processes them in parallel.

**Inputs (per sample):**

| Input | Description |
| --- | --- |
| `gvcf` | The per-sample GVCF file to reblock. |
| `gvcf_index` | The GVCF's `.tbi` index. |
| `ref_fasta`, `ref_fasta_index`, `ref_dict` | hg38 reference files (pre-configured in this workspace). |
| `cloud_provider` | Set to `"gcp"` when running in Terra on Google Cloud. |

**Output:**

| Output | Description |
| --- | --- |
| `reblocked_gvcf` | The reblocked GVCF, named `<sample>.rb.g.vcf.gz`. |
| `reblocked_gvcf_index` | Its `.tbi` index. |

**To run it in Terra:**

1. Go to the **Workflows** tab and select the **ReblockGVCF** workflow.
2. Configure it to run on your **sample** data table so each row (sample) is reblocked.
3. Point `gvcf` / `gvcf_index` at your per-sample GVCF columns; leave the reference inputs at their pre-configured values; set `cloud_provider` to `"gcp"`.
4. Launch. When it finishes, the reblocked GVCF paths are written back to your data table — use these in the next step.

---

## Step 1 — Build a Sample Map

The JointGenotyping workflow does not take a list of GVCFs directly. Instead it takes a **sample map** (`sample_name_map`): a **tab-separated file with two columns and no header**:

```
Sample_A    gs://your-bucket/Sample_A.rb.g.vcf.gz
Sample_B    gs://your-bucket/Sample_B.rb.g.vcf.gz
Sample_C    gs://your-bucket/Sample_C.rb.g.vcf.gz
```

- **Column 1:** the sample name (must be unique across the cohort).
- **Column 2:** the cloud (`gs://`) path to that sample's GVCF (use the **reblocked** GVCF if you ran Step 0).

This lets the pipeline know the sample names without downloading every GVCF header. You can create this file with a helper workflow (e.g. a **Generate-Sample-Map** workflow that reads sample names and GVCF paths from a data table) or by hand for small cohorts. Upload the resulting `sample_map.tsv` to the workspace bucket and record its path in your `sample_set` data table.

---

## Step 2 — Run JointGenotyping

**What the pipeline does (high level):** it imports all your GVCFs into a joint database, genotypes them together, filters the resulting variants, collects quality-control metrics, and (optionally) checks sample fingerprints to guard against sample swaps. The result is a **single filtered multi-sample VCF** in which every sample has a genotype at every retained site; sites that fail filtering are **not removed** but are flagged in the VCF `FILTER` column.

For a full, step-by-step description of every task and tool the pipeline runs, see the **[JointGenotyping pipeline documentation](https://broadinstitute.github.io/warp/docs/Pipelines/JointGenotyping_Pipeline/README)**. The essentials for running it are below.

### Key inputs

| Input | Description |
| --- | --- |
| `sample_name_map` | The tab-separated sample map from [Step 1](#step-1--build-a-sample-map). |
| `callset_name` | A name for your cohort; used to name the output files (e.g. `Kinshasa_cohort_2026`). |
| `unpadded_intervals_file` | Genome intervals over which to call (pre-configured for hg38). |
| Reference & resource files | hg38 FASTA, dbSNP, HapMap, Omni, 1000G, Mills, Axiom, `haplotype_database`, `eval_interval_list` — **all pre-configured in this workspace**. |
| `small_disk` / `medium_disk` / `large_disk` / `huge_disk` | Disk sizes; scale with cohort size. Use the values in the example configuration as a starting point. |

The workflow accepts **nested inputs** (`allowNestedInputs: true`), so you generally only need to set `sample_name_map` and `callset_name` and leave everything else at the workspace defaults.

### Variant filtering: VQSR (default)

By default the pipeline filters variants using **VQSR** (Variant Quality Score Recalibration), the long-standing [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612) method. VQSR builds a statistical model of what "real" SNPs and indels look like using trusted resource files (HapMap, Omni, 1000G, Mills, etc.), then flags likely-artifact sites in the `FILTER` column. **For most users, leave the defaults as-is** and VQSR will be used.

### Advanced filtering options (optional)

You can change the filtering/genotyping strategy if you have a specific reason to:

- **VETS** (Variant Extract-Train-Score) — a newer filtering approach. Enable with `run_vets = true`. Requires a `targets_interval_list`.
- **GnarlyGenotyper** — a faster, "quick and dirty" joint genotyping method. Enable with `use_gnarly_genotyper = true`.

If you are new to joint calling, we recommend starting with the **default VQSR** path and only exploring these once you are comfortable.

### To run it in Terra

1. Go to the **Workflows** tab and select the **JointGenotyping** workflow.
2. Choose your `sample_set` data table as the root entity and select your cohort.
3. Confirm `sample_name_map` and `callset_name` are set; leave the reference inputs and filtering options at their defaults for a standard VQSR run.
4. Launch the workflow and monitor progress in the **Job History** tab.

---

## Step 3 — Outputs

The pipeline produces the following outputs (written to the workspace bucket):

| Output | Filename | Description |
| --- | --- | --- |
| `output_vcfs` | `<callset_name>.vcf.gz` (small cohorts) or `<callset_name>.filtered.<idx>.vcf.gz` (large cohorts) | The filtered, joint-called multi-sample VCF(s). **This is what you load into seqr.** |
| `output_vcf_indices` | `<callset_name>.vcf.gz.tbi` (and/or per-shard `.tbi`) | Index files for the VCF(s). |
| `detail_metrics_file` | `<callset_name>.variant_calling_detail_metrics` | Per-sample variant-calling QC metrics (Picard). |
| `summary_metrics_file` | `<callset_name>.variant_calling_summary_metrics` | Cohort-level variant-calling QC metrics (Picard). |
| `output_intervals` | interval list | The intervals used for this run. |
| `crosscheck_fingerprint_check` | `<callset_name>.fingerprintcheck` | (Optional) fingerprint check to detect sample swaps. |

> **Important for large cohorts:** for callsets with **more than ~1000 samples**, the pipeline leaves the VCF **split into multiple shards** (`output_vcfs` is an array of per-interval files) rather than one combined file. seqr expects a **single** joint VCF, so for large cohorts you will need to **merge the shards** into one VCF before loading. For typical cohorts (tens to hundreds of samples), the pipeline already produces a single `<callset_name>.vcf.gz` and no merge is needed.

---

## Step 4 — Handing off to seqr

[**seqr**](https://seqr.broadinstitute.org/) is an open-source platform for rare-disease genomics that lets you search and interpret variants across a family or cohort. seqr takes a **joint-called (multi-sample) VCF** as its input and runs its **own annotation pipeline** on it — so the filtered VCF from Step 3 is the correct starting point.

**Things to know before loading into seqr:**

- **One joint VCF.** seqr expects a single multi-sample VCF. For small/medium cohorts this pipeline already gives you one; for very large cohorts, merge the shards first (see the note above).
- **Reference build.** seqr supports **GRCh38 (hg38)**, which is exactly what this pipeline produces.
- **Sample IDs.** The sample names in your VCF (which come from your sample map) should match the sample IDs in the pedigree/metadata you provide to seqr.
- **Annotation is handled by seqr.** You do **not** need to pre-annotate the VCF (e.g. with VEP) — seqr's loading pipeline does this.

**How to load into seqr** depends on how you access it:

- **Via [AnVIL](https://anvilproject.org/):** you can request that a joint-called VCF be loaded into seqr. See the [seqr on AnVIL video tutorial](https://www.youtube.com/watch?v=hnAI55cqGbk).
- **As a Broad [CMG](https://cmg.broadinstitute.org/) / [GREGoR](https://gregorconsortium.org/) collaborator.**
- **On a self-hosted seqr instance** (seqr is [open source on GitHub](https://github.com/broadinstitute/seqr)).

For general information, training videos, and contact details, see the [seqr website](https://seqr.broadinstitute.org/) and the [seqr publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9903206). *(This section is intentionally general; once the specific seqr loading path for this project is confirmed, we can add exact, step-by-step handoff instructions here.)*

---

## Time and Cost Estimates

<!-- TODO: fill in with measured Terra benchmarks for the example cohort. -->

Cost and runtime depend heavily on cohort size and whether GVCFs were reblocked. Benchmarks below are placeholders to be filled in once measured.

### ReblockGVCF (per sample)

| Sample | Time | Cost $ |
| --- | --- | --- |
| _Whole genome GVCF (example)_ | 1h 27m | $0.04 |

### JointGenotyping

| Cohort | # Samples | Reblocked? | Time | Cost $ |
| --- | --- | --- | --- | --- |
| WGS_example_data | 66 | Yes | 6h 6m | $20.29 |
| WGS_example_data | 66 | No | 8h 57m | $93.02 |

For guidance on controlling cloud costs in Terra, see [this article](https://support.terra.bio/hc/en-us/articles/360029748111).

---

## Additional Resources

- **JointGenotyping full pipeline documentation:** [broadinstitute.github.io/warp/docs/Pipelines/JointGenotyping_Pipeline/README](https://broadinstitute.github.io/warp/docs/Pipelines/JointGenotyping_Pipeline/README)
- **JointGenotyping WDL source & changelog:** [WARP repository](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/dna_seq/germline/joint_genotyping)
- **ReblockGVCF WDL source:** [reblocking folder](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/dna_seq/germline/joint_genotyping/reblocking)
- **Reblocking background:** [GATK reblocking article](https://gatk.broadinstitute.org/hc/en-us/articles/4414594365467) and the [WARP reblocking blog post](https://broadinstitute.github.io/warp/blog/Nov21_ReblockedGVCF)
- **GATK Best Practices — Germline short variant discovery:** [GATK documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932)
- **seqr:** [seqr.broadinstitute.org](https://seqr.broadinstitute.org/)
- **Terra support:** [support.terra.bio](https://support.terra.bio/hc/en-us)

---

## Contact Information

- For workspace questions and feedback, email the Broad pipelines team at warp-pipelines-help@broadinstitute.org, or [file an issue in WARP](https://github.com/broadinstitute/warp/issues).
- For GATK tool questions, see the [GATK forum](https://gatk.broadinstitute.org/hc/en-us/community/topics).
- For seqr questions, see the [seqr website](https://seqr.broadinstitute.org/) or email seqr@broadinstitute.org.

---

## Citing the JointGenotyping Pipeline

If you use the JointGenotyping pipeline in your research, please consider citing:

> Degatano, K., Awdeh, A., Cox III, R.S., Dingman, W., Grant, G., Khajouei, F., Kiernan, E., Konwar, K., Mathews, K.L., Palis, K., et al. Warp Analysis Research Pipelines: Cloud-optimized workflows for biological data processing and reproducible analysis. *Bioinformatics* 2025; btaf494. https://doi.org/10.1093/bioinformatics/btaf494

---

## License

**Copyright Broad Institute, 2026 | BSD-3**
All code referenced in this workspace is released under the BSD-3 license (full text at https://github.com/broadinstitute/warp/blob/develop/LICENSE). Programs called by the workflows may be subject to different licenses; users are responsible for confirming they are authorized to run all tools.

---

## Workspace Change Log

| Date | Change | Author |
| --- | --- | --- |
| _TBD_ | Initial Joint Calling tutorial for the Kinshasa consulting engagement. | _TBD_ |
