version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow FilterAFComparisonTable {
	input {
		File table
	}

	call FilterAFComparisonTable_t {
		input:
			table = table
	}

	output {
		File filtered_table = FilterAFComparisonTable_t.filtered_table
	}
}

task FilterAFComparisonTable_t {
	input {
		File table
	}

	command <<<
		Rscript -<<"EOF"
			library(dplyr)
			library(readr)

			thousand_genomes_afr_haplotypes <- 2 * 662
			thousand_genomes_amr_haplotypes <- 2 * 347
			thousand_genomes_eas_haplotypes <- 2 * 504
			thousand_genomes_eur_haplotypes <- 2 * 503
			thousand_genomes_haplotypes <- thousand_genomes_afr_haplotypes + thousand_genomes_amr_haplotypes + thousand_genomes_eas_haplotypes + thousand_genomes_eur_haplotypes

			filter_sites <- function(x, pos) {
				ret <- x %>% rowwise() %>% filter(!is.na(gnomad.controls_AF_afr)) %>% mutate(thousand_g_n = AFR_AF * thousand_genomes_afr_haplotypes
																						+ AMR_AF * thousand_genomes_amr_haplotypes
																						+ EAS_AF * thousand_genomes_eas_haplotypes
																						+ EUR_AF * thousand_genomes_eur_haplotypes,
																						gnomad_p = (gnomad.controls_AF_afr * thousand_genomes_afr_haplotypes
																									+ gnomad.controls_AF_amr * thousand_genomes_amr_haplotypes
																									+ gnomad.controls_AF_eas * thousand_genomes_eas_haplotypes
																									+ gnomad.controls_AF_nfe * thousand_genomes_eur_haplotypes)/
																							(thousand_genomes_afr_haplotypes
																							+ thousand_genomes_amr_haplotypes
																							+ thousand_genomes_eas_haplotypes
																							+ thousand_genomes_eur_haplotypes
																							),
																						thousand_g_n_wo_eur = AFR_AF * thousand_genomes_afr_haplotypes
																						+ AMR_AF * thousand_genomes_amr_haplotypes
																						+ EAS_AF * thousand_genomes_eas_haplotypes,
																						gnomad_p_wo_eur = (gnomad.controls_AF_afr * thousand_genomes_afr_haplotypes
																											+ gnomad.controls_AF_amr * thousand_genomes_amr_haplotypes
																											+ gnomad.controls_AF_eas * thousand_genomes_eas_haplotypes)/
																							thousand_genomes_haplotypes)
				ret <- ret %>% filter(thousand_g_n > 0) %>% filter(gnomad_p > 0) %>% rowwise() %>% mutate(p_value = as.numeric(binom.test(round(thousand_g_n), thousand_genomes_haplotypes, gnomad_p)$p.value),
																										p_value_wo_eur = as.numeric(binom.test(round(thousand_g_n_wo_eur), thousand_genomes_haplotypes - thousand_genomes_eur_haplotypes, gnomad_p_wo_eur)$p.value))

				ret <- ret %>% filter(p_value < 1e-10 & p_value_wo_eur < 1e-10)
				ret
			}

			filtered_sites <- read_tsv_chunked("~{table}", DataFrameCallback$new(filter_sites), col_types = cols(CHROM=col_character()))

			write_tsv(filtered_sites, "filteredAFComparisonSites.tsv")
		EOF
	>>>

	runtime {
		docker: "rocker/tidyverse"
		disks: "local-disk 100 HDD"
		memory: "32 GB"
	}

	output {
		File filtered_table = "filteredAFComparisonSites.tsv"
	}
}