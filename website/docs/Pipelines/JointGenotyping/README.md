---
sidebar_position: 1
slug: /Pipelines/JointGenotyping_Pipeline/README
---

# JointGenotyping Overview

| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [JointGenotyping_v1.6.10](https://github.com/broadinstitute/warp/releases) | February, 2024 | Elizabeth Kiernan & Kaylee Mathews | Please [file an issue in WARP](https://github.com/broadinstitute/warp/issues). |

## Introduction to the JointGenotyping workflow

The [JointGenotyping workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/dna_seq/germline/joint_genotyping/JointGenotyping.wdl) is an open-source, cloud-optimized pipeline that implements joint variant calling, filtering, and (optional) fingerprinting.

The pipeline can be configured to run using one of the following GATK joint genotyping methods:

* **[GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/21905118377755)** (default method) performs joint genotyping on GVCF files stored in GenomicsDB and pre-called with HaplotypeCaller.
* **[GnarlyGenotyper](https://gatk.broadinstitute.org/hc/en-us/articles/21904951112091)** performs scalable, “quick and dirty” joint genotyping on a set of GVCF files stored in GenomicsDB and pre-called with HaplotypeCaller.

The pipeline can be configured to run using one of the following GATK variant filtering techniques:

* **[Variant Quality Score Recalibration (VQSR)](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612)** (default method) uses the VariantRecalibrator and ApplyVQSR tools to filter variants according to [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932).
* **Variant Extract-Train-Score (VETS)** uses the ExtractVariantAnnotations, TrainVariantAnnotationsModel, and ScoreVariantAnnotations tools called in the [VETS subworkflow](https://github.com/broadinstitute/gatk/blob/master/scripts/vcf_site_level_filtering_wdl/JointVcfFiltering.wdl) to score variant annotations.

The pipeline takes in a sample map file listing GVCF files produced by HaplotypeCaller in GVCF mode and produces a filtered VCF file (with index) containing genotypes for all samples present in the input VCF files. All sites that are present in the input VCF file are retained. Filtered sites are annotated as such in the FILTER field. If you are new to VCF files, see the [file type specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

The JointGenotyping pipeline can be adapted to run on Microsoft Azure instead of Google Cloud. For more information, see the [azure-warp-joint-calling GitHub repository](https://github.com/broadinstitute/azure-warp-joint-calling).

## Set-up

### JointGenotyping Installation and Requirements

To download the latest JointGenotyping release, see the release tags prefixed with "JointGenotyping" on the WARP [releases page](https://github.com/broadinstitute/warp/releases). All JointGenotyping pipeline releases are documented in the [JointGenotyping changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/dna_seq/germline/joint_genotyping/JointGenotyping.changelog.md). 

To search releases of this and other pipelines, use the WARP command-line tool [Wreleaser](https://github.com/broadinstitute/warp/tree/master/wreleaser).

If you’re running a JointGenotyping workflow version prior to the latest release, the accompanying documentation for that release may be downloaded with the source code on the WARP [releases page](https://github.com/broadinstitute/warp/releases) (see the folder `website/docs/Pipelines/JointGenotyping`).

The JointGenotyping pipeline can be deployed using [Cromwell](https://cromwell.readthedocs.io/en/stable/), a GA4GH-compliant, flexible workflow management system that supports multiple computing platforms. The workflow can also be run in [Terra](https://app.terra.bio), a cloud-based analysis platform. The Terra [Whole-Genome-Analysis-Pipeline](https://app.terra.bio/#workspaces/warp-pipelines/Whole-Genome-Analysis-Pipeline) and [Exome-Analysis-Pipeline](https://app.terra.bio/#workspaces/warp-pipelines/Exome-Analysis-Pipeline) workspaces contain the JointGenotyping pipeline, as well as workflows for preprocessing, initial variant calling, and sample map generation, workflow configurations, required reference data and other inputs, and example testing data.

### Inputs

The JointGenotyping workflow inputs are specified in JSON configuration files. Example configuration files can be found in the [test_inputs](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/dna_seq/germline/joint_genotyping/test_inputs) folder in the WARP repository.

#### Default joint calling input descriptions

The table below describes the pipeline inputs that apply when the pipeline is run with default parameters and uses GenotypeGVCFs for joint calling and VQSR for variant filtering:

| Parameter name | Description | Type |
| --- | --- | --- |
| unpadded_intervals_file | Describes the intervals for which VCF output will be written; exome data will have different captures/targets. | File |
| callset_name | Identifier for the group of VCF files used for joint calling. | String |
| sample_name_map | Path to file containing the sample names and the cloud location of the individual GVCF files. | String | 
| ref_fasta | Reference FASTA file used for joint calling; must agree with reference for `unpadded_intervals_file`. | File |
| ref_fasta_index | Index for reference FASTA file used for joint calling; must agree with reference for `unpadded_intervals_file`. | File |
| ref_dict | Reference dictionary file used for joint calling; must agree with reference for `unpadded_intervals_file`. | File |
| dbsnp_vcf | Resource VCF file containing common SNPs and indels used for annotating the VCF file after joint calling. | File |
| dbsnp_vcf_index | Index for `dbsnp_vcf`. | File |
| snp_recalibration_tranche_values | Set of sensitivity levels used when running the pipeline using VQSR; value should match estimated sensitivity of truth resource passed as `hapmap_resource_vcf` to the [SNPsVariantRecalibratorCreateModel](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) and [SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) tasks; filter cutoff based on sensitivity to common variants (more sensitivity = more false positives); required when `run_vets` is “false”. | Array[String] |
| snp_recalibration_annotation_values | Features used for filtering model (annotations in VCF file); all allele-specific versions. | Array[String] |
| indel_recalibration_tranche_values | Set of sensitivity levels used when running the pipeline using VQSR; value should match estimated sensitivity of truth resource passed as `mills_resource_vcf` to the [IndelsVariantRecalibrator](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task; filter cutoff based on sensitivity to common variants (more sensitivity = more false positives); required when `run_vets` is “false”. | Array[String] |
| indel_recalibration_annotation_values | Features used for filtering model when running the pipeline using VQSR; required when `run_vets` is “false”. | Array[String] |
| eval_interval_list | Subset of the unpadded intervals file used for metrics.  | File |
| hapmap_resource_vcf | Used for SNP variant recalibration; see the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) for more information. | File |
| hapmap_resource_vcf_index | Used for SNP variant recalibration; see the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) for more information. | File |
| omni_resource_vcf | Used for SNP recalibration; see the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) for more information. | File |
| omni_resource_vcf_index | Used for SNP recalibration; see the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) for more information. | File | 
| one_thousand_genomes_resource_vcf | Used for SNP recalibration; see the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) for more information. | File |
| one_thousand_genomes_resource_vcf_index | Used for SNP recalibration; see the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) for more information. | File |
| mills_resource_vcf | Used for indel variant recalibration; see the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) for more information. | File |
| mills_resource_vcf_index | Used for indel variant recalibration; see the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) for more information. | File |
| axiomPoly_resource_vcf | Used for indel variant recalibration; see the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) for more information. | File |
| axiomPoly_resource_vcf_index | Used for indel variant recalibration; see the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) for more information. | File |
| dbsnp_resource_vcf | Optional file used for SNP/indel variant recalibration; set to `dbsnp_vcf` by default; see the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) for more information. | File |
| dbsnp_resource_vcf_index | Optional file used for SNP/indel variant recalibration; set to `dbsnp_vcf_index` by default; see the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811) for more information. | File |
| excess_het_threshold | Optional float used for hard filtering joint calls; phred-scaled p-value; set to `54.69` by default to cut off quality scores greater than a z-score of -4.5 (p-value of 3.4e-06). | Float |
| vqsr_snp_filter_level | Used for applying the recalibration model when running the pipeline using VQSR; required when `run_vets` is “false”. | Float |
| vqsr_indel_filter_level | Used for applying the recalibration model when running the pipeline using VQSR; required when `run_vets` is “false”. | Float |
| snp_vqsr_downsampleFactor | The downsample factor used for SNP variant recalibration if the number of GVCF files is greater than the ` snps_variant_recalibration_threshold` when running the pipeline using VQSR; required when `run_vets` is “false”. | Int |
| top_level_scatter_count | Optional integer used to determine how many files the input interval list should be split into; default will split the interval list into 2 files. | Int |
| gather_vcfs | Optional boolean; “true” is used for small callsets containing less than 100,000 GVCF files. | Boolean |
| snps_variant_recalibration_threshold | Optional integer that sets the threshold for the number of callset VCF files used to perform recalibration on a single file; if the number of VCF files exceeds the threshold, variants will be downsampled to enable parallelization; default is “500000”. | Int | 
| rename_gvcf_samples | Optional boolean describing whether GVCF samples should be renamed; default is “true”. | Boolean |
| unbounded_scatter_count_scale_factor | Optional float used to scale the scatter count when `top_level_scatter_count` is not provided as input; default is “0.15”. | Float |
| use_allele_specific_annotations | Optional boolean used for SNP and indel variant recalibration when running the pipeline using VQSR; set to “true” by default. | Boolean |


#### GnarlyGenotyper joint calling input descriptions

The table below describes the additional pipeline inputs that apply when the pipeline is run with GnarlyGenotyper for joint calling:

| Parameter name | Description | Type |
| --- | --- | --- |
| gnarly_scatter_count | Optional integer used to determine how many files to split the interval list into when using GnarlyGenotyper; default is “10”. | Int |
| use_gnarly_genotyper | Optional boolean describing whether GnarlyGenotyper should be used; default is “false”. | Boolean |


#### VETS variant filtering input descriptions

The table below describes the additional pipeline inputs that apply when the pipeline is run with VETS for variant filtering:

| Parameter name | Description | Type |
| --- | --- | --- |
| targets_interval_list | Describes the intervals for which the filtering model will be trained when running the pipeline using VETS; for more details, see the associated [README](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/bge_exome_calling_regions.v1.1.interval_list.README.md); required when `run_vets` is “true”. | File |
| run_vets | Optional boolean used to describe whether the pipeline will use VQSR (`run_vets = false`) or VETS (`run_vets = true`) to create the filtering model; default is “false”. | Boolean |


#### Fingerprinting input descriptions

The table below describes the pipeline inputs that apply to fingerprinting:

| Parameter name | Description | Type |
| --- | --- | --- |
| haplotype_database | Haplotype reference used for fingerprinting (see the CrosscheckFingerprints task). | File |
| cross_check_fingerprints | Optional boolean describing whether or not the pipeline should check fingerprints; default is “true”. | Boolean |
| scatter_cross_check_fingerprints | Optional boolean describing whether `CrossCheckFingerprintsScattered` or `CrossCheckFingerprintsSolo` should be run; default is “false” and `CrossCheckFingerprintsSolo` will be run. | Boolean |

#### Runtime parameter input descriptions

The table below describes the pipeline inputs used for setting runtime parameters of tasks:

| Parameter name | Description | Type |
| --- | --- | --- |
| small_disk | Disk size; dependent on cohort size; requires user input; see example JSON configuration files found in the WARP [test_inputs](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/dna_seq/germline/joint_genotyping/test_inputs) folder for recommendations. | Int |
| medium_disk | Disk size; dependent on cohort size; requires user input; see example JSON configuration files found in the WARP [test_inputs](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/dna_seq/germline/joint_genotyping/test_inputs) folder for recommendations. | Int |
| large_disk | Disk size; dependent on cohort size; requires user input; see example JSON configuration files found in the WARP [test_inputs](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/dna_seq/germline/joint_genotyping/test_inputs) folder for recommendations. | Int | 
| huge_disk | Disk size; dependent on cohort size; requires user input; see example JSON configuration files found in the WARP [test_inputs](https://github.com/broadinstitute/warp/tree/master/pipelines/wdl/dna_seq/germline/joint_genotyping/test_inputs) folder for recommendations. | Int |


## JointGenotyping tasks and tools

The [JointGenotyping workflow](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/dna_seq/germline/joint_genotyping/JointGenotyping.wdl) imports individual "tasks," also written in WDL script, from the WARP [tasks folder](https://github.com/broadinstitute/warp/tree/master/tasks/broad). 

Overall, the JointGenotyping workflow:

1. Splits the input interval list and imports GVCF files.
1. Performs joint genotyping using GATK GenotypeGVCFs (default) or GnarlyGenotyper.
1. Creates single site-specific VCF and index files.
1. Creates and applies a variant filtering model using GATK VQSR (default) or VETS.
1. Collects variant calling metrics.
1. Checks fingerprints (optional).

The tasks and tools used in the JointGenotyping workflow are detailed in the table below. 

To see specific tool parameters, select the task WDL link in the table; then find the task and view the `command {}` section of the task in the WDL script. To view or use the exact tool software, see the task's Docker image which is specified in the task WDL `# runtime values` section as `String docker =`.

| Task | Tool | Software | Description | 
| --- | --- | --- | --- | 
| [CheckSamplesUnique](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | bash | bash | Checks that there are more than 50 unique samples in `sample_name_map`. | 
| [SplitIntervalList](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | SplitIntervals | [GATK](https://gatk.broadinstitute.org/hc/en-us) | Splits the unpadded interval list for scattering. |
| [ImportGVCFs](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | GenomicsDBImport | [GATK](https://gatk.broadinstitute.org/hc/en-us) | Imports single-sample GVCF files into GenomicsDB before joint genotyping. |
| [SplitIntervalList as GnarlyIntervalScatterDude](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | SplitIntervals | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `use_gnarly_genotyper` is “true” (default is “false”), splits the unpadded interval list for scattering; otherwise, this task is skipped. |
| [GnarlyGenotyper](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | GnarlyGenotyper | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `use_gnarly_genotyper` is “true” (default is “false”), performs scalable, “quick and dirty” joint genotyping on a set of GVCF files stored in GenomicsDB; otherwise, this task is skipped. |
| [GatherVcfs as TotallyRadicalGatherVcfs](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | GatherVcfsCloud | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `use_gnarly_genotyper` is “true” (default is “false”), compiles the site-specific VCF files generated for each interval into one VCF output and index; otherwise, this task is skipped. | 
| [GenotypeGVCFs](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | GenotypeGVCFs | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `use_gnarly_genotyper` is “false” (default is “false”), performs joint genotyping on GVCF files stored in GenomicsDB; otherwise this task is skipped. |
| [HardFilterAndMakeSitesOnlyVcf](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | VariantFiltration, MakeSitesOnlyVcf | [GATK](https://gatk.broadinstitute.org/hc/en-us) | Uses the VCF files to hard filter the variant calls; outputs a VCF file with the site-specific (but not genotype) information. | 
| [GatherVcfs as SitesOnlyGatherVcf](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | GatherVcfsCloud | [GATK](https://gatk.broadinstitute.org/hc/en-us) | Compiles the site-specific VCF files generated for each interval into one VCF output file and index. | 
| [JointVcfFiltering as TrainAndApplyVETS](https://github.com/broadinstitute/gatk/blob/master/scripts/vcf_site_level_filtering_wdl/JointVcfFiltering.wdl) | ExtractVariantAnnotations, TrainVariantAnnotationsModel, ScoreVariantAnnotations | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `run_vets` is “true” (default is “false”), calls the `JointVcfFiltering.wdl` subworkflow to extract variant-level annotations, trains a model for variant scoring, and scores variants; otherwise, this task is skipped. | 
| [IndelsVariantRecalibrator](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | VariantRecalibrator | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `run_vets` is “false” (default is “false”), uses the compiled VCF file to build a recalibration model to score indel variant quality; produces a recalibration table. | 
| [SNPsVariantRecalibratorCreateModel](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | VariantRecalibrator | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `run_vets` is “false” (default is “false”) and the number of input GVCF files is greater than `snps_variant_recalibration_threshold`, builds a recalibration model to score variant quality; otherwise this task is skipped. | 
| [SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | VariantRecalibrator | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `run_vets` is “false” (default is “false”) and the number of input GVCF files is greater than `snps_variant_recalibration_threshold`, builds a scattered recalibration model to score variant quality; otherwise this task is skipped. | 
| [Tasks.GatherTranches as SNPGatherTranches](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | GatherTranches | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `run_vets` is “false” (default is “false”) and the number of input GVCF files is greater than `snps_variant_recalibration_threshold`, gathers tranches into a single file; otherwise this task is skipped. | 
| [SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | VariantRecalibrator | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `run_vets` is “false” (default is “false”) and the number of input GVCF files is not greater than `snps_variant_recalibration_threshold`, builds a recalibration model to score variant quality; otherwise this task is skipped. | 
| [ApplyRecalibration](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | ApplyVQSR | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `run_vets` is “false” (default is “false”), scatters the site-specific VCF file and applies a filtering threshold. |
| [CollectVariantCallingMetrics as CollectMetricsSharded](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | CollectVariantCallingMetrics | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If the callset has at least 1000 GVCF files, returns detail and summary metrics for each of the scattered VCF files. If the number is small, will return metrics for a merged VCF file produced in the `GatherVcfs as FinalGatherVcf` task (listed below). | 
| [GatherVcfs as FinalGatherVcf](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | GatherVcfsCloud | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If the callset has fewer than 1000 GVCF files, compiles the VCF files prior to collecting metrics in the `CollectVariantCallingMetrics as CollectMetricsOnFullVcf` task (listed below). |
| [CollectVariantCallingMetrics as CollectMetricsOnFullVcf](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | CollectVariantCallingMetrics | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If the callset has fewer than 1000 GVCF files, returns metrics for the merged VCF file produced in the `GatherVcfs as FinalGatherVcf` task. |
| [GatherVariantCallingMetrics](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | AccumulateVariantCallingMetrics | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If the callset has at least 1000 GVCF files, gathers metrics produced for each VCF file. |
| [GetFingerprintingIntervalIndices](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | IntervalListTools | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `cross_check_fingerprints` is “true” (default is “true”) and `scatter_cross_check_fingerprints` is “true” (default is “false”), gets and sorts indices for fingerprint intervals; otherwise the task is skipped. |
| [GatherVcfs as GatherFingerprintingVcfs](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | GatherVcfsCloud | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `cross_check_fingerprints` is “true” (default is “true”) and `scatter_cross_check_fingerprints` is “true” (default is “false”), compiles the fingerprint VCF files; otherwise the task is skipped. |
| [SelectFingerprintSiteVariants](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | SelectVariants | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `cross_check_fingerprints` is “true”  (default is “true”)and `scatter_cross_check_fingerprints` is “true” (default is “false”), selects variants from the fingerprint VCF file; otherwise the task is skipped. |
| [PartitionSampleNameMap](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | bash | bash | If `cross_check_fingerprints` is “true” (default is “true”) and `scatter_cross_check_fingerprints` is “true” (default is “false”), partitions the sample name map and files are scattered by the partition; otherwise the task is skipped. |
| [CrossCheckFingerprint as CrossCheckFingerprintsScattered](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | CrosscheckFingerprints | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `cross_check_fingerprints` is “true” (default is “true”) and `scatter_cross_check_fingerprints` is “true” (default is “false”), checks fingerprints for the VCFs in the scattered partitions and produces a metrics file; otherwise the task is skipped. |
| [GatherPicardMetrics as GatherFingerprintingMetrics](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | bash | bash | If `cross_check_fingerprints` is “true” (default is “true”) and `scatter_cross_check_fingerprints` is “true” (default is “false”), combines the fingerprint metrics files into a single metrics file; otherwise the task is skipped. |
| [CrossCheckFingerprint as CrossCheckFingerprintSolo](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) | CrosscheckFingerprints | [GATK](https://gatk.broadinstitute.org/hc/en-us) | If `cross_check_fingerprints` is “true” (default is “true”) and `scatter_cross_check_fingerprints` is “false” (default is “false”), checks fingerprints for the single VCF file and produces a metrics file; otherwise the task is skipped. |

#### 1. Splits the input interval list and imports GVCF files

The [SplitIntervalList](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task uses GATK’s SplitIntervals tool to split the input interval list into two or more interval files. The number of output interval files can be specified using the `top_level_scatter_count` input parameter or by specifying `unbounded_scatter_count_scale_factor`, which will scale the number of output files based on the number of input GVCF files.

The [ImportGVCFs](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task uses GATK’s GenomicsDBImport tool and the input sample map file to import single-sample GVCF files into GenomicsDB before joint genotyping.

#### 2. Performs joint genotyping using GATK GenotypeGVCFs (default) or GnarlyGenotyper

**GenotypeGVCFs (default)**

When `use_gnarly_genotyper` is “false”, the [GenotypeGVCFs](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task uses GATK’s GenotypeGVCFs tool to perform joint genotyping on GVCF files stored in GenomicsDB that have been pre-called with HaplotypeCaller.

**GnarlyGenotyper**

When `use_gnarly_genotyper` is “true”, the [SplitIntervalList as GnarlyIntervalScatterDude](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task splits the unpadded interval list for scattering using GATK’s SplitIntervals tool. The output is used as input for the [GnarlyGenotyper](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task which performs joint genotyping on the set of GVCF files and outputs an array of VCF and index files using the GnarlyGenotyper tool. Those VCF and index files are gathered in the next task, [GatherVcfs as TotallyRadicalGatherVcfs](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl), which uses the GatherVcfsCloud tool.

#### 3. Creates single site-specific VCF and index files

The [HardFilterAndMakeSitesOnlyVcf](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task takes in the output VCF and index files produced by either GnarlyGenotyper or GenotypeGVCFs. The task uses the `excess_het_threshold` input value to hard filter the variant calls using GATK’s VariantFiltration tool. After filtering, the site-specific VCF files are generated from the filtered VCF files by removing all sample-specific genotype information, leaving only the site-level summary information at each site.

Next, the site-specific VCF and index files for each interval are gathered into a single site-specific VCF and index file by the [GatherVcfs as SitesOnlyGatherVcf](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task, which uses the GatherVcfsCloud tool.

#### 4. Creates and applies a variant filtering model using GATK VQSR (default) or VETS

**VQSR (default)**

If `run_vets` is “false”, the [IndelsVariantRecalibrator](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task takes in the site-specific VCF and index files generated in [Step 3](#3-creates-single-site-specific-vcf-and-index-files) and uses GATK’s VariantRecalibrator tool to perform the first step of the Variant Quality Score Recalibration (VQSR) technique of filtering variants. The tool builds a model to be used to score and filter indels and produces a recalibration table as output.

After building the indel filtering model, the workflow uses the VariantRecalibrator tool to build a model to be used to score and filter SNPs. If the number of input GVCF files is greater than `snps_variant_recalibration_threshold`, the [SNPsVariantRecalibratorCreateModel](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl), [SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl), and [Tasks.GatherTranches as SNPGatherTranches](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) tasks are called to scatter the site-specific VCF and index files, build the SNP model, and gather scattered tranches into a single file. If the number of input GVCF files is less than `snps_variant_recalibration_threshold`, the [SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task is called to build the SNP model.

The [ApplyRecalibration](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task uses GATK’s ApplyVQSR tool to scatter the site-specific VCF file, apply the indel and SNP filtering models, and output a recalibrated VCF and index file. 

**VETS**

If `run_vets` is “true”, the [JointVcfFiltering as TrainAndApplyVETS](https://github.com/broadinstitute/gatk/blob/master/scripts/vcf_site_level_filtering_wdl/JointVcfFiltering.wdl) task takes in the hard filtered and site-specific VCF and index files generated in [Step 3](#3-creates-single-site-specific-vcf-and-index-files) and calls the `JointVcfFiltering.wdl` subworkflow. This workflow uses the Variant Extract-Train-Score (VETS) algorithm to extract variant-level annotations, train a filtering model, and score variants based on the model. The subworkflow uses the GATK ExtractVariantAnnotations, TrainVariantAnnotationsModel, and ScoreVariantAnnotations tools to create extracted and scored VCF and index files. The output VCF and index files are not filtered by the score assigned by the model. The score is included in the output VCF files in the INFO field as an annotation called “SCORE”.

The VETS algorithm trains the model only over target regions, rather than including exon tails which can lead to poor-quality data. However, the model is applied everywhere including the exon tails.

#### 5. Collects variant calling metrics

Summary and per-sample metrics are collected using Picard’s CollectVariantCallingMetrics tool. For large callsets (at least 1000 GVCF files), the workflow calls the [CollectVariantCallingMetrics as CollectMetricsSharded](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) followed by the [GatherVariantCallingMetrics](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task to compute and gather the variant calling metrics into single output files. For small callsets (less than 1000 GVCF files), the workflow calls the [GatherVcfs as FinalGatherVcf](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task followed by the [CollectVariantCallingMetrics as CollectMetricsOnFullVcf](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task to first compile the VCF files and then compute the variant calling metrics. Detail and summary metrics files are produced as outputs of these tasks.

#### 6. Checks fingerprints (optional)

If `cross_check_fingerprints` is “true”, the workflow will use Picard to determine the likelihood that the input and output data were generated from the same individual to verify that the pipeline didn’t swap any of the samples during processing. The [SelectFingerprintSiteVariants](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/JointGenotypingTasks.wdl) task uses GATK’s SelectVariants tool to select variants in the site-specific VCF file based on the variants present in the `haplotype_database` and outputs a fingerprint VCF and index file. Next, the workflow cross-checks the fingerprints and creates an output metrics file using the CrosscheckFingerprints tool.

## Outputs

The following table lists the output variables and files produced by the pipeline.

| Output name | Filename, if applicable | Output format and description |
| ------ | ------ | ------ |
| detail_metrics_file | `<callset_name>.variant_calling_detail_metrics` | Detail metrics file produced using Picard. |
| summary_metrics_file | `<callset_name>.variant_calling_summary_metrics` | Summary metrics file produced using Picard. |
| output_vcfs | `<callset_name>.vcf.gz` or `<callset_name>.filtered.<idx>.vcf.gz` | Array of all site-specific output VCF files. |
| output_vcf_indices | `<callset_name>.vcf.gz.tbi` or `<callset_name>.filtered.<idx>.vcf.gz.tbi` | Array of all output VCF index files. |
| output_intervals | `scatterDir/<output_intervals_files>` | Interval list file produced by the workflow. |
| crosscheck_fingerprint_check | `<callset_name>.fingerprintcheck` | Fingerprint metrics | Optional output file containing fingerprint metrics. |

## Versioning and testing

All JointGenotyping pipeline releases are documented in the [JointGenotyping changelog](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/dna_seq/germline/joint_genotyping/JointGenotyping.changelog.md) and tested using [plumbing and scientific test data](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/dna_seq/germline/joint_genotyping/test_data_overview.md). To learn more about WARP pipeline testing, see [Testing Pipelines](https://broadinstitute.github.io/warp/docs/About_WARP/TestingPipelines).

## Citing the JointGenotyping Pipeline

If you use the JointGenotyping Pipeline in your research, please consider citing our preprint:

Degatano, K., Awdeh, A., Cox III, R.S., Dingman, W., Grant, G., Khajouei, F., Kiernan, E., Konwar, K., Mathews, K.L., Palis, K., et al. Warp Analysis Research Pipelines: Cloud-optimized workflows for biological data processing and reproducible analysis. Bioinformatics 2025; btaf494. https://doi.org/10.1093/bioinformatics/btaf494

## Feedback

Please help us make our tools better by [filing an issue in WARP](https://github.com/broadinstitute/warp/issues); we welcome pipeline-related suggestions or questions.