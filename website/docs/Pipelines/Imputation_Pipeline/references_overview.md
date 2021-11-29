# Creating the Imputation Pipeline's Modified 1000 Genomes Reference

## Background
Initial tests of the Imputation workflow followed by assessments of polygenic risk score revealed that disease risk scores were lower when computed from imputed array data as opposed to whole-genome sequencing data (see figures below), possibly because of incorrectly genotyped sites in the 1000 Genomes (1000G) reference panel. 

![](imputed_vs_wgs_scores_original_1kg-1.png)

![](imputed_vs_wgs_scores_original_1kg-2.png)

This systematic difference was found to be due to a single site (10:104952499) which has a relatively high weight in the CAD weights file and looks to have been incorrectly genotyped in 1000 Genomes. This site has an allele fraction of 0.72 in 1000 Genomes, but only 0.086 in the [gnomAD V2 reference](https://gnomad.broadinstitute.org/). 


As a result, the 1000G reference files were modified for the Imputation pipeline as described below. You can view the original, unmodified 1000G VCFs [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). 

## Reference modification 
To remove putative incorrect sites from the 1000G reference panel, allele frequencies were compared between it and the gnomAD V2 reference panel. First, the [BuildAFComparisonTable workflow](https://github.com/broadinstitute/warp/tree/develop/scripts/BuildAFComparisonTable.wdl) was used to create a table of the allele frequencies for both reference panels. Then, the [FilterAFComparisonTable workflow](https://github.com/broadinstitute/warp/tree/develop/scripts/FilterAFComparisonTable.wdl) was applied to compare the observed number of alleles in 1000G to the expected number of alleles set by the gnomAD V2 reference using a two-sided binomial p-value. If both p-values were less than 1e-10, then the site was flagged as incorrect. After identifying the putative incorrect sites, the [RemoveBadSitesById workflow](https://github.com/broadinstitute/warp/tree/develop/scripts/RemoveBadSitesById.wdl) was used to remove them, generating a cleaned 1000G reference panel. 

This cleaning removes 359,369, or about 0.8% of sites from 1000 Genomes reference. In the histogram below, only sites that were were flagged as incorrect are shown. The vast majority of flagged sites have p-values that are much lower than the threshold.

![](method_2_p-value_histogram-1.png)

As can be seen below, using this "cleaned" 1000G removes the systematic difference between imputed and WGS scores.

![](imputed_vs_wgs_scores_cleaned_method_2_1kg-1.png)
![](imputed_vs_wgs_scores_cleaned_method_2_1kg-2.png)

The comparison below shows that the improvement is due only to the removal of sites that were poorly imputed with the original 1000G reference, the remaining sites are not being imputed at a higher quality than they were with the original 1000G reference. 

![](cleaned_vs_original_compare_to_gnomad_af_method_correlations-1.png)

A public copy of the cleaned reference can be found at gs://broad-gotc-test-storage/imputation/1000G_reference_panel/ as shown in the Imputation workflow's [example configuration file](https://github.com/broadinstitute/warp/blob/master/pipelines/broad/arrays/imputation/example_inputs.json) (JSON).

## Questions 
For questions or additional information about the Imputation pipeline's reference generation, email [Chris Kachullis](mailto:ckachuli@broadinstitute.org).

