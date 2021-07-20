version 1.0

import "../../../../structs/imputation/ImputationStructs.wdl" as structs

workflow ImputationPipeline {
	input {
		Int chunkLength = 25000000
		Int chunkOverlaps = 5000000 # this is the padding that will be added to the beginning and end of each chunk to reduce edge effects

		# You can either input a multisample VCF or an array of single sample VCFs
		# The pipeline will just merge the single sample VCFs into one multisample VCF
		# and then impute the multisample VCF
		# If you want to run a single sample VCF, set the multi_sample_vcf input to the
		# single sample VCF
		File? multi_sample_vcf
		File? multi_sample_vcf_index
		Array[File]? single_sample_vcfs
		Array[File]? single_sample_vcf_indices

		Boolean perform_extra_qc_steps # these are optional additional extra QC steps from Amit's group that should only be
		# run for large sample sets, especially a diverse set of samples (it's further limiting called at sites to 95% and by HWE)
		File ref_dict # for reheadering / adding contig lengths in the header of the ouptut VCF, and calculating contig lengths
		Array[ReferencePanelContig] referencePanelContigs
		String genetic_maps_eagle = "/genetic_map_hg19_withX.txt.gz" # this is for Eagle, it is in the docker image
		String output_callset_name # the output callset name
		Boolean split_output_to_single_sample = false
		File haplotype_database
		Int merge_ssvcf_mem = 3 # the memory allocation for MergeSingleSampleVcfs (in GiB)
	}
	# Docker images here
	String bcftools_docker_tag = "us.gcr.io/broad-dsde-methods/imputation_bcftools_vcftools_docker:v1.0.0" #"us.gcr.io/broad-dsde-methods/bcftools:v1.9-1-deb_cv1"
	String bcftools_vcftools_docker_tag = "us.gcr.io/broad-dsde-methods/imputation_bcftools_vcftools_docker:v1.0.0" #"us.gcr.io/broad-dsde-methods/skwalker-imputation:with_vcftools" # TODO: use a public one (not suse bcftools biocontainers also has vcftools )
	String gatk_docker_tag = "us.gcr.io/broad-gatk/gatk:4.1.9.0"
	String minimac4_docker_tag = "us.gcr.io/broad-dsde-methods/imputation-minimac-docker:v1.0.0" #"us.gcr.io/broad-dsde-methods/skwalker-imputation:test" # this has the exact version of minimac and eagle we want to perfectly match Michigan Server
	String eagle_docker_tag = "us.gcr.io/broad-dsde-methods/imputation_eagle_docker:v1.0.0"
	String ubuntu_docker_tag = "ubuntu:20.04"
	String rtidyverse_docker_tag = "rocker/tidyverse:4.1.0"

	if (defined(single_sample_vcfs)) {
		call MergeSingleSampleVcfs {
			input:
				input_vcfs = select_first([single_sample_vcfs, []]),
				input_vcf_indices = select_first([single_sample_vcf_indices, []]),
				output_vcf_basename = "merged_input_samples",
				bcftools_docker = bcftools_docker_tag,
				mem = merge_ssvcf_mem
		}
	}

	File vcf_to_impute = select_first([multi_sample_vcf, MergeSingleSampleVcfs.output_vcf])
	File vcf_index_to_impute = select_first([multi_sample_vcf_index, MergeSingleSampleVcfs.output_vcf_index])

	call SetIDs {
		input:
			vcf = vcf_to_impute,
			vcfIndex = vcf_index_to_impute,
			output_basename = "input_samples_with_variant_ids",
			bcftools_docker = bcftools_docker_tag
	}

	call SortIds as SortIdsVcfToImpute {
		input:
			vcf = SetIDs.output_vcf,
			basename = "input_samples_with_variant_ids_sorted",
			bcftools_docker = bcftools_docker_tag
	}

	call ExtractIDs as ExtractIdsVcfToImpute {
		input:
			vcf = SortIdsVcfToImpute.output_vcf,
			output_basename = "imputed_sites",
			bcftools_docker = bcftools_docker_tag
	}

	call CountSamples {
		input:
			vcf = vcf_to_impute,
			bcftools_docker = bcftools_docker_tag
	}

	scatter (referencePanelContig in referencePanelContigs) {
		call CalculateChromsomeLength {
			input:
				ref_dict = ref_dict,
				chrom = referencePanelContig.contig,
				gatk_docker = gatk_docker_tag
		}

		Float chunkLengthFloat = chunkLength
		Int num_chunks = ceil(CalculateChromsomeLength.chrom_length / chunkLengthFloat)

		scatter (i in range(num_chunks)) {
			String chunk_contig = referencePanelContig.contig
			Int start = (i * chunkLength) + 1
			Int startWithOverlaps = if (start - chunkOverlaps < 1) then 1 else start - chunkOverlaps
			Int end = if (CalculateChromsomeLength.chrom_length < ((i + 1) * chunkLength)) then CalculateChromsomeLength.chrom_length else ((i + 1) * chunkLength)
			Int endWithOverlaps = if (CalculateChromsomeLength.chrom_length < end + chunkOverlaps) then CalculateChromsomeLength.chrom_length else end + chunkOverlaps

			call GenerateChunk {
				input:
					vcf = vcf_to_impute,
					vcf_index = vcf_index_to_impute,
					start = startWithOverlaps,
					end = endWithOverlaps,
					chrom = referencePanelContig.contig,
					basename = "chrom_" + referencePanelContig.contig + "_chunk_" + i,
					gatk_docker = gatk_docker_tag
			}

			if (perform_extra_qc_steps) {
				call OptionalQCSites {
					input:
						input_vcf = GenerateChunk.output_vcf,
						input_vcf_index = GenerateChunk.output_vcf_index,
						output_vcf_basename =  "chrom_" + referencePanelContig.contig + "_chunk_" + i,
						bcftools_vcftools_docker = bcftools_vcftools_docker_tag
				}
			}

			call CountChunks {
				input:
					vcf = select_first([OptionalQCSites.output_vcf,  GenerateChunk.output_vcf]),
					vcf_index = select_first([OptionalQCSites.output_vcf_index, GenerateChunk.output_vcf_index]),
					panel_vcf = referencePanelContig.vcf,
					panel_vcf_index = referencePanelContig.vcf_index,
					gatk_docker = gatk_docker_tag
			}
			call CheckChunks {
				input:
					vcf = select_first([OptionalQCSites.output_vcf,  GenerateChunk.output_vcf]),
					vcf_index = select_first([OptionalQCSites.output_vcf_index, GenerateChunk.output_vcf_index]),
					panel_vcf = referencePanelContig.vcf,
					panel_vcf_index = referencePanelContig.vcf_index,
					var_in_original = CountChunks.var_in_original,
					var_in_reference = CountChunks.var_in_reference,
					bcftools_docker = bcftools_docker_tag
			}

			if (CheckChunks.valid) {

				call PrePhaseVariantsEagle {
					input:
						dataset_bcf = CheckChunks.valid_chunk_bcf,
						dataset_bcf_index = CheckChunks.valid_chunk_bcf_index,
						reference_panel_bcf = referencePanelContig.bcf,
						reference_panel_bcf_index = referencePanelContig.bcf_index,
						chrom = referencePanelContig.contig,
						genetic_map_file = genetic_maps_eagle,
						eagle_docker = eagle_docker_tag,
						start = startWithOverlaps,
						end = endWithOverlaps
				}

				call minimac4 {
					input:
						ref_panel = referencePanelContig.m3vcf,
						phased_vcf = PrePhaseVariantsEagle.dataset_prephased_vcf,
						prefix = "chrom" + "_chunk_" + i +"_imputed",
						chrom = referencePanelContig.contig,
						minimac4_docker = minimac4_docker_tag,
						start = start,
						end = end,
						window = chunkOverlaps
				}

				call AggregateImputationQCMetrics {
					input:
						infoFile = minimac4.info,
						nSamples = CountSamples.nSamples,
						basename = output_callset_name + "chrom_" + referencePanelContig.contig + "_chunk_" + i,
						rtidyverse_docker = rtidyverse_docker_tag
				}

				call UpdateHeader {
					input:
						vcf = minimac4.vcf,
						vcf_index = minimac4.vcf_index,
						ref_dict = ref_dict,
						basename = "chrom" + "_chunk_" + i +"_imputed",
						gatk_docker = gatk_docker_tag
				}

				call SeparateMultiallelics {
					input:
						original_vcf = UpdateHeader.output_vcf,
						original_vcf_index = UpdateHeader.output_vcf_index,
						output_basename = "chrom" + referencePanelContig.contig + "_chunk_" + i +"_imputed",
						bcftools_docker = bcftools_docker_tag
				}

				call RemoveSymbolicAlleles {
					input:
						original_vcf = SeparateMultiallelics.output_vcf,
						original_vcf_index = SeparateMultiallelics.output_vcf_index,
						output_basename = "chrom" + referencePanelContig.contig + "_chunk_" + i +"_imputed",
						gatk_docker = gatk_docker_tag
				}

				call SortIds {
					input:
						vcf = RemoveSymbolicAlleles.output_vcf,
						basename = "chrom" + referencePanelContig.contig + "_chunk_" + i +"_imputed",
						bcftools_docker = bcftools_docker_tag
				}
			}
		}
		Array[File] aggregatedImputationMetrics = select_all(AggregateImputationQCMetrics.aggregated_metrics)
		Array[File] chromosome_vcfs = select_all(SortIds.output_vcf)
		Array[File] chromosome_vcf_indices = select_all(SortIds.output_vcf_index)
	}

	Array[File] phased_vcfs = flatten(chromosome_vcfs)
	Array[File] phased_vcf_indices = flatten(chromosome_vcf_indices)

	call GatherVcfs {
		input:
			input_vcfs = phased_vcfs,
			input_vcf_indices = phased_vcf_indices,
			output_vcf_basename = output_callset_name,
			gatk_docker = gatk_docker_tag
	}

	call ExtractIDs {
		input:
			vcf = GatherVcfs.output_vcf,
			output_basename = "imputed_sites",
			bcftools_docker = bcftools_docker_tag
	}

	call FindSitesFileTwoOnly {
		input:
			file1 = ExtractIDs.ids,
			file2 = ExtractIdsVcfToImpute.ids,
			ubuntu_docker = ubuntu_docker_tag
	}

	call SelectVariantsByIds {
		input:
			vcf = SortIdsVcfToImpute.output_vcf,
			ids = FindSitesFileTwoOnly.missing_sites,
			basename = "imputed_sites_to_recover",
			gatk_docker = gatk_docker_tag
	}

	call RemoveAnnotations {
		input:
			vcf = SelectVariantsByIds.output_vcf,
			basename = "imputed_sites_to_recover_annotations_removed",
			bcftools_docker = bcftools_docker_tag
	}

	call InterleaveVariants {
		input:
			vcfs = [RemoveAnnotations.output_vcf, GatherVcfs.output_vcf],
			basename = output_callset_name,
			gatk_docker = gatk_docker_tag
	}

	call MergeImputationQCMetrics {
		input:
			metrics = flatten(aggregatedImputationMetrics),
			basename = output_callset_name,
			rtidyverse_docker = rtidyverse_docker_tag
	}


	call StoreChunksInfo {
		input:
			chroms = flatten(chunk_contig),
			starts = flatten(start),
			ends = flatten(end),
			vars_in_array = flatten(CountChunks.var_in_original),
			vars_in_panel = flatten(CountChunks.var_in_reference),
			valids = flatten(CheckChunks.valid),
			basename = output_callset_name,
			rtidyverse_docker = rtidyverse_docker_tag
	}

	call CrosscheckFingerprints {
		input:
			firstInputs = if (defined(multi_sample_vcf)) then select_all([multi_sample_vcf]) else select_first([single_sample_vcfs]),
			firstInputIndices = if (defined(multi_sample_vcf)) then select_all([multi_sample_vcf_index]) else select_first([single_sample_vcf_indices]),
			secondInputs = [InterleaveVariants.output_vcf],
			secondInputIndices = [InterleaveVariants.output_vcf_index],
			haplotypeDatabase = haplotype_database,
			basename = output_callset_name,
			gatk_docker = gatk_docker_tag
	}

	if (split_output_to_single_sample) {
		call SplitMultiSampleVcf {
			input:
				multiSampleVcf = InterleaveVariants.output_vcf,
				bcftools_docker = bcftools_docker_tag
		}

		call CrosscheckFingerprints as CrosscheckFingerprintsSplit {
			input:
				firstInputs = if (defined(multi_sample_vcf)) then select_all([multi_sample_vcf]) else select_first([single_sample_vcfs]),
				firstInputIndices = if (defined(multi_sample_vcf)) then select_all([multi_sample_vcf_index]) else select_first([single_sample_vcf_indices]),
				secondInputs = SplitMultiSampleVcf.single_sample_vcfs,
				secondInputIndices = SplitMultiSampleVcf.single_sample_vcf_indices,
				haplotypeDatabase = haplotype_database,
				basename = output_callset_name + ".split",
				gatk_docker = gatk_docker_tag
		}
	}


	output {
		Array[File]? imputed_single_sample_vcfs = SplitMultiSampleVcf.single_sample_vcfs
		Array[File]? imputed_single_sample_vcf_indices = SplitMultiSampleVcf.single_sample_vcf_indices
		File imputed_multisample_vcf = InterleaveVariants.output_vcf
		File imputed_multisample_vcf_index = InterleaveVariants.output_vcf_index
		File aggregated_imputation_metrics = MergeImputationQCMetrics.aggregated_metrics
		File chunks_info = StoreChunksInfo.chunks_info
		File failed_chunks = StoreChunksInfo.failed_chunks
		Int n_failed_chunks = StoreChunksInfo.n_failed_chunks
		File crosscheck = CrosscheckFingerprints.crosscheck
		File? crosscheck_split = CrosscheckFingerprintsSplit.crosscheck
	}
}

task RemoveDuplicates {
	input {
		File input_vcfs
		File input_vcf_indices
		String output_basename
		String bcftools_docker
	}
	Int disk_size = ceil(size([input_vcfs, input_vcf_indices], "GB"))
	command <<<
		bcftools norm -d both ~{input_vcfs} | bgzip -c > ~{output_basename}.vcf.gz
		bcftools index -t ~{output_basename}.vcf.gz
	>>>
	runtime {
		docker: bcftools_docker
		memory: "3 GiB"
		disks: "local-disk " + disk_size + " HDD"
	}
	output {
		File output_vcf = "~{output_basename}.vcf.gz"
		File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
	}
}

task CalculateChromsomeLength {
	input {
		File ref_dict
		Int chrom
		String gatk_docker
	}
	command {
		grep -P "SN:~{chrom}\t" ~{ref_dict} | sed 's/.*LN://' | sed 's/\t.*//'
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk 100 HDD"
		memory: "2 GB"
	}
	output {
		Int chrom_length = read_int(stdout())
	}
}

task GenerateChunk {
	input {
		Int start
		Int end
		String chrom
		String basename
		String vcf
		String vcf_index
		Int disk_size = 400 # not sure how big the disk size needs to be since we aren't downloading the entire VCF here
		String gatk_docker
	}
	command {
		gatk SelectVariants -V ~{vcf} --select-type-to-include SNP --max-nocall-fraction 0.1 -xl-select-type SYMBOLIC \
		--select-type-to-exclude MIXED --restrict-alleles-to BIALLELIC -L ~{chrom}:~{start}-~{end} -O ~{basename}.vcf.gz --exclude-filtered true
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: "8 GB"
	}
	parameter_meta {
		vcf: {
				 description: "vcf",
				 localization_optional: true
			 } vcf_index: {
							  description: "vcf index",
							  localization_optional: true
						  }
	}
	output {
		File output_vcf = "~{basename}.vcf.gz"
		File output_vcf_index = "~{basename}.vcf.gz.tbi"
	}
}

task CountChunks {
	input {
		File vcf
		File vcf_index
		File panel_vcf
		File panel_vcf_index
		Int disk_size =ceil(2*size([vcf, vcf_index, panel_vcf, panel_vcf_index], "GB"))
		String gatk_docker
	}
	command <<<
		echo $(gatk CountVariants -V ~{vcf}  | sed 's/Tool returned://') > var_in_original
		echo $(gatk  CountVariants -V ~{vcf} -L ~{panel_vcf}  | sed 's/Tool returned://') > var_in_reference
		echo ${var_in_reference} " * 2 - " ${var_in_original} "should be greater than 0 AND " ${var_in_reference} "should be greater than 3"
	>>>
	output {
		Int var_in_original = read_int("var_in_original")
		Int var_in_reference = read_int("var_in_reference")
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: "4 GB"
	}

}
task CheckChunks {
	input {
		File vcf
		File vcf_index
		File panel_vcf
		File panel_vcf_index
		Int var_in_original
		Int var_in_reference
		Int disk_size =ceil(2*size([vcf, vcf_index, panel_vcf, panel_vcf_index], "GB"))
		String bcftools_docker
	}
	command <<<
		if [ $(( ~{var_in_reference} * 2 - ~{var_in_original})) -gt 0 ] && [ ~{var_in_reference} -gt 3 ]; then
		echo true > valid_file.txt
		bcftools convert -Ob ~{vcf} > valid_variants.bcf
		bcftools index -f valid_variants.bcf
		else
		echo false > valid_file.txt
		fi
	>>>
	output {
		File? valid_chunk_bcf = "valid_variants.bcf"
		File? valid_chunk_bcf_index = "valid_variants.bcf.csi"
		Boolean valid = read_boolean("valid_file.txt")
	}
	runtime {
		docker: bcftools_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: "4 GB"
	}

}

task PrePhaseVariantsEagle {
	input {
		File? dataset_bcf
		File? dataset_bcf_index
		File reference_panel_bcf
		File reference_panel_bcf_index
		String chrom
		String genetic_map_file
		Int start
		Int end
		String eagle_docker
	}
	Int disk_size = ceil(3 * size([dataset_bcf, reference_panel_bcf, dataset_bcf_index, reference_panel_bcf_index], "GB"))
	command <<<
		/eagle  \
		--vcfTarget ~{dataset_bcf}  \
		--vcfRef ~{reference_panel_bcf} \
		--geneticMapFile ~{genetic_map_file} \
		--outPrefix pre_phased_~{chrom} \
		--vcfOutFormat z \
		--bpStart ~{start} --bpEnd ~{end} --allowRefAltSwap
	>>>
	output {
		File dataset_prephased_vcf="pre_phased_~{chrom}.vcf.gz"
	}
	runtime {
		docker: eagle_docker
		memory: "32 GB"
		cpu: "8"
		disks: "local-disk " + disk_size + " HDD"
	}
}

task minimac4 {
	input {
		File ref_panel
		File phased_vcf
		String prefix
		String chrom
		Int start
		Int end
		String minimac4_docker
		Int window
	}
	command <<<
		/Minimac4 --refHaps ~{ref_panel} --haps ~{phased_vcf} --start ~{start} --end ~{end} --window ~{window} \
		--chr ~{chrom} --noPhoneHome --format GT,DS,GP --allTypedSites --prefix ~{prefix} --minRatio 0.00001
		if [ ! -f ~{prefix}.dose.vcf.gz.tbi ]
		then
		bcftools index -t ~{prefix}.dose.vcf.gz
		fi
	>>>
	output {
		File vcf = "~{prefix}.dose.vcf.gz"
		File vcf_index = "~{prefix}.dose.vcf.gz.tbi"
		File info = "~{prefix}.info"
	}
	runtime {
		docker: minimac4_docker
		memory: "4 GB"
		cpu: "1"
		disks: "local-disk 100 HDD"
	}
}

task MergeVCFs {
	input {
		Array[File] input_vcfs # these are all valid
		String output_vcf_basename
		String bcftools_docker
	}

	Int disk_size = ceil(3*size(input_vcfs, "GB"))

	command <<<
		bcftools concat ~{sep=' ' input_vcfs} -Oz -o ~{output_vcf_basename}.vcf.gz
		bcftools index -t ~{output_vcf_basename}.vcf.gz
	>>>
	runtime {
		docker: bcftools_docker
		memory: "3 GiB"
		disks: "local-disk " + disk_size + " HDD"
	}
	output {
		File output_vcf = "~{output_vcf_basename}.vcf.gz"
		File output_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
	}
}

task GatherVcfs {
	input {
		Array[File] input_vcfs # these are all valid
		Array[File] input_vcf_indices # these are all valid
		String output_vcf_basename
		String gatk_docker
	}

	Int disk_size = ceil(3*size(input_vcfs, "GB"))

	command <<<
		gatk GatherVcfs -I ~{sep=' -I ' input_vcfs} -O ~{output_vcf_basename}.vcf.gz # I don't think this creates an output index file
		if [ ! -f ~{output_vcf_basename}.vcf.gz.tbi ]
		then
		gatk IndexFeatureFile -I ~{output_vcf_basename}.vcf.gz
		fi
	>>>
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: "16 GB"
	}
	output {
		File output_vcf = "~{output_vcf_basename}.vcf.gz"
		File output_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
	}
}

task UpdateHeader {
	input {
		File vcf
		File vcf_index
		File ref_dict
		String basename
		Int disk_size = ceil(4*(size(vcf, "GB") + size(vcf_index, "GB"))) + 20
		String gatk_docker
	}
	command <<<

		## update the header of the merged vcf
		gatk UpdateVCFSequenceDictionary --source-dictionary ~{ref_dict} --output ~{basename}.vcf.gz \
		--replace -V ~{vcf} --disable-sequence-dictionary-validation
	>>>
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: "8 GiB"
	}
	output {
		File output_vcf = "~{basename}.vcf.gz"
		File output_vcf_index = "~{basename}.vcf.gz.tbi"
	}
}

task RemoveSymbolicAlleles {
	input {
		File original_vcf
		File original_vcf_index
		String output_basename
		Int disk_size = ceil(3*(size(original_vcf, "GB") + size(original_vcf_index, "GB")))
		String gatk_docker
	}
	command {
		gatk SelectVariants -V ~{original_vcf} -xl-select-type SYMBOLIC -O ~{output_basename}.vcf.gz
	}
	output {
		File output_vcf = "~{output_basename}.vcf.gz"
		File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: "4 GB"
	}
}

task SeparateMultiallelics {
	input {
		File original_vcf
		File original_vcf_index
		String output_basename
		Int disk_size =  ceil(2*(size(original_vcf, "GB") + size(original_vcf_index, "GB")))
		String bcftools_docker
	}
	command {
		bcftools norm -m - ~{original_vcf} -Oz -o ~{output_basename}.vcf.gz
		bcftools index -t ~{output_basename}.vcf.gz
	}
	output {
		File output_vcf = "~{output_basename}.vcf.gz"
		File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
	}
	runtime {
		docker: bcftools_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: "4 GB"
	}
}

task SplitSample {
	input {
		File vcf
		File vcf_index
		String sample
		Int disk_size = ceil(2*(size(vcf, "GB") + size(vcf_index, "GB")))
		String gatk_docker
	}
	command {
		gatk SelectVariants -V ~{vcf} -sn ~{sample} -O ~{sample}.vcf.gz
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: "9 GB"
	}
	output {
		File output_gzipped_vcf = "~{sample}.vcf.gz"
		File output_gzipped_vcf_index = "~{sample}.vcf.gz.tbi"
	}
}

task OptionalQCSites {
	input {
		File input_vcf
		File input_vcf_index
		String output_vcf_basename
		String bcftools_vcftools_docker
	}
	Int disk_size = ceil(2*(size(input_vcf, "GB") + size(input_vcf_index, "GB")))

	command <<<
		# site missing rate < 5% ; hwe p > 1e-6
		vcftools --gzvcf ~{input_vcf}  --max-missing 0.05 --hwe 0.000001 --recode -c | bgzip -c > ~{output_vcf_basename}.vcf.gz
		bcftools index -t ~{output_vcf_basename}.vcf.gz
	>>> # Note: this is necessary because vcftools doesn't have a way to output a zipped vcf, nor a way to index one (hence needing to use bcf).

	runtime {
		docker: bcftools_vcftools_docker
		memory: "16 GiB"
		disks: "local-disk " + disk_size + " HDD"
	}
	output {
		File output_vcf = "~{output_vcf_basename}.vcf.gz"
		File output_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
	}
}

task SortIds {
	input {
		File vcf
		String basename
		Int disk_space =  3*ceil(size(vcf, "GB"))
		String bcftools_docker
	}

	command <<<
		# what better way is there to do this I really don't know
		zcat ~{vcf} | awk -v OFS='\t' '{split($3, n, ":"); if ( !($1 ~ /^"#"/) && n[4] < n[3])  $3=n[1]":"n[2]":"n[4]":"n[3]; print $0}' | bgzip -c > ~{basename}.vcf.gz
		bcftools index -t ~{basename}.vcf.gz

	>>>

	output {
		File output_vcf = "~{basename}.vcf.gz"
		File output_vcf_index = "~{basename}.vcf.gz.tbi"

	}

	runtime {
		docker: bcftools_docker
		disks: "local-disk " + disk_space + " HDD"
		memory: "16 GB"
	}
}

task MergeSingleSampleVcfs {
	input {
		Array[File] input_vcfs
		Array[File] input_vcf_indices
		String output_vcf_basename
		String bcftools_docker
		Int mem
	}

	Int disk_size = 3* ceil(size(input_vcfs, "GB") + size(input_vcf_indices, "GB")) + 20

	command <<<
		bcftools merge ~{sep=' ' input_vcfs} -Oz -o ~{output_vcf_basename}.vcf.gz
		bcftools index -t ~{output_vcf_basename}.vcf.gz
	>>>

	runtime {
		docker: bcftools_docker
		memory: mem + " GiB"
		disks: "local-disk " + disk_size + " HDD"
	}
	output {
		File output_vcf = "~{output_vcf_basename}.vcf.gz"
		File output_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
	}
}

task CountSamples {
	input {
		File vcf
		String bcftools_docker
	}

	Int disk_size = 100 + ceil(size(vcf, "GB"))

	command <<<
		bcftools query -l ~{vcf} | wc -l >> nSamples.txt
	>>>

	runtime {
		docker: bcftools_docker
		memory: "3 GiB"
		disks: "local-disk " + disk_size + " HDD"
	}

	output {
		Int nSamples = read_int("nSamples.txt")
	}
}

task AggregateImputationQCMetrics {
	input {
		File infoFile
		Int nSamples
		String basename
		String rtidyverse_docker
	}

	Int disk_size = 100 + ceil(size(infoFile, "GB"))

	command <<<
		Rscript -<< "EOF"
		library(dplyr)
		library(readr)
		library(purrr)
		library(ggplot2)

		sites_info <- read_tsv("~{infoFile}")

		nSites <- sites_info %>% nrow()
		nSites_with_var <- sites_info %>% filter(MAF >= 0.3/(2*~{nSamples} - 0.7)) %>% nrow()
		nSites_high_r2 <- sites_info %>% filter(Rsq>0.3) %>% nrow()

		aggregated_metrics <- tibble(total_sites=nSites, total_sites_with_var=nSites_with_var, total_sites_r2_gt_0.3=nSites_high_r2,)

		write_tsv(aggregated_metrics, "~{basename}_aggregated_imputation_metrics.tsv")

		EOF
	>>>

	runtime {
		docker: rtidyverse_docker
		disks : "local-disk " + disk_size + " HDD"
		preemptible : 3
	}

	output {
		File aggregated_metrics = "~{basename}_aggregated_imputation_metrics.tsv"
	}
}

task StoreChunksInfo {
	input {
		Array[String] chroms
		Array[Int] starts
		Array[Int] ends
		Array[Int] vars_in_array
		Array[Int] vars_in_panel
		Array[Boolean] valids
		String basename
		String rtidyverse_docker
	}

	command <<<
		Rscript -<< "EOF"
		library(dplyr)
		library(readr)

		chunk_info <- tibble(chrom = c("~{sep='", "' chroms}"), start = c("~{sep='", "' starts}"), ends = c("~{sep='", "' ends}"), vars_in_array = c("~{sep='", "' vars_in_array}"), vars_in_panel = c("~{sep='", "' vars_in_panel}"), chunk_was_imputed = as.logical(c("~{sep='", "' valids}")))
		failed_chunks <- chunk_info %>% filter(!chunk_was_imputed) %>% select(-chunk_was_imputed)
		n_failed_chunks <- nrow(failed_chunks)
		write_tsv(chunk_info, "~{basename}_chunk_info.tsv")
		write_tsv(failed_chunks, "~{basename}_failed_chunks.tsv")
		write(n_failed_chunks, "n_failed_chunks.txt")
		EOF
	>>>

	runtime {
		docker: rtidyverse_docker
		preemptible : 3
	}

	output {
		File chunks_info = "~{basename}_chunk_info.tsv"
		File failed_chunks = "~{basename}_failed_chunks.tsv"
		Int n_failed_chunks = read_int("n_failed_chunks.txt")
	}
}

task MergeImputationQCMetrics {
	input {
		Array[File] metrics
		String basename
		String rtidyverse_docker
	}

	Int disk_size = 100 + ceil(size(metrics, "GB"))

	command <<<
		Rscript -<< "EOF"
		library(dplyr)
		library(readr)
		library(purrr)
		library(ggplot2)

		metrics <- list("~{sep='", "' metrics}") %>% map(read_tsv) %>% reduce(`+`) %>% mutate(frac_sites_r2_gt_0.3=total_sites_r2_gt_0.3/total_sites, frac_sites_with_var_r2_gt_0.3=total_sites_r2_gt_0.3/total_sites_with_var)

		write_tsv(metrics, "~{basename}_aggregated_imputation_metrics.tsv")

		EOF
	>>>

	runtime {
		docker: rtidyverse_docker
		disks : "local-disk " + disk_size + " HDD"
		preemptible : 3
	}

	output {
		File aggregated_metrics = "~{basename}_aggregated_imputation_metrics.tsv"
	}
}

task SetIDs {
	input {
		File vcf
		File vcfIndex
		String output_basename
		String bcftools_docker
	}

	Int disk_size = 100 + ceil(2.2 * size(vcf, "GB"))

	command <<<
		bcftools annotate ~{vcf} --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' -Oz -o ~{output_basename}.vcf.gz
		bcftools index -t ~{output_basename}.vcf.gz
	>>>

	runtime {
		docker: bcftools_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: "4 GB"
	}

	output {
		File output_vcf = "~{output_basename}.vcf.gz"
		File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
	}
}

task ExtractIDs {
	input {
		File vcf
		String output_basename
		Int disk_size = 2*ceil(size(vcf, "GB")) + 100
		String bcftools_docker
	}

	command <<<
		bcftools query -f "%ID\n" ~{vcf} -o ~{output_basename}.ids
	>>>
	output {
		File ids = "~{output_basename}.ids"
	}
	runtime {
		docker: bcftools_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: "4 GB"
	}
}

task SelectVariantsByIds {
	input {
		File vcf
		File ids
		String basename
		String gatk_docker
	}

	Int disk_size = ceil(1.2*size(vcf, "GB")) + 100

	parameter_meta {
		vcf: {
				 description: "vcf",
				 localization_optional: true
			 }
	}

	command <<<
		cp ~{ids} sites.list
		gatk SelectVariants -V ~{vcf} --exclude-filtered --keep-ids sites.list -O ~{basename}.vcf.gz
	>>>

	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " SSD"
		memory: "16 GB"
	}

	output {
		File output_vcf = "~{basename}.vcf.gz"
		File output_vcf_index = "~{basename}.vcf.gz.tbi"
	}
}

task RemoveAnnotations {
	input {
		File vcf
		String basename
		String bcftools_docker
	}

	Int disk_size = ceil(2.2*size(vcf, "GB")) + 100

	command <<<
		bcftools annotate ~{vcf} -x FORMAT,INFO -Oz -o ~{basename}.vcf.gz
		bcftools index -t ~{basename}.vcf.gz
	>>>

	runtime {
		docker: bcftools_docker
		memory: "3 GiB"
		disks: "local-disk " + disk_size + " HDD"
	}

	output {
		File output_vcf = "~{basename}.vcf.gz"
		File output_vcf_index = "~{basename}.vcf.gz.tbi"
	}
}

task InterleaveVariants {
	input {
		Array[File] vcfs
		String basename
		String gatk_docker
	}

	Int disk_size = ceil(3.2*size(vcfs, "GB")) + 100

	command <<<
		gatk MergeVcfs -I ~{sep=" -I " vcfs} -O ~{basename}.vcf.gz
	>>>


	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " SSD"
		memory: "16 GB"
	}

	output {
		File output_vcf = "~{basename}.vcf.gz"
		File output_vcf_index = "~{basename}.vcf.gz.tbi"
	}
}

task FindSitesFileTwoOnly {
	input {
		File file1
		File file2
		String ubuntu_docker
	}

	Int disk_size = ceil(size(file1, "GB") + 2*size(file2, "GB")) + 100

	command <<<
		comm -13 <(sort ~{file1} | uniq) <(sort ~{file2} | uniq) > missing_sites.ids
	>>>

	runtime {
		docker: ubuntu_docker
		disks: "local-disk " + disk_size + " SSD"
		memory: "16 GB"
	}

	output {
		File missing_sites = "missing_sites.ids"
	}
}

task SplitMultiSampleVcf {
	input {
		File multiSampleVcf
		Int mem = 8
		Int bcftools_docker
	}

	Int disk_size = ceil(3*size(multiSampleVcf, "GB")) + 100

	command <<<
		mkdir out_dir
		bcftools +split ~{multiSampleVcf} -Oz -o out_dir
		for vcf in out_dir/*.vcf.gz; do
		bcftools index -t $vcf
		done
	>>>

	runtime {
		docker: bcftools_docker
		disks: "local-disk " + disk_size + " SSD"
		memory: mem + " GB"
	}

	output {
		Array[File] single_sample_vcfs = glob("out_dir/*.vcf.gz")
		Array[File] single_sample_vcf_indices = glob("out_dir/*.vcf.gz.tbi")
	}
}

task CrosscheckFingerprints {
	input {
		Array[File] firstInputs
		Array[File] secondInputs
		Array[File] firstInputIndices
		Array[File] secondInputIndices
		File haplotypeDatabase
		String basename
		Int mem = 8
		String gatk_docker
	}

	Int disk_size = ceil(1.2*(size(firstInputs, "GB") + size(secondInputs, "GB") + size(haplotypeDatabase, "GB"))) + 100

	command <<<
		# add links to ensure correctly located indices
		array_vcfs=( ~{sep=" " firstInputs} )
		array_indices=( ~{sep=" " firstInputIndices} )
		for i in ${!array_vcfs[@]}; do
		ln -s ${array_indices[i]} $(dirname ${array_vcfs[i]})
		done

		array_vcfs2=( ~{sep=" " secondInputs} )
		array_indices2=( ~{sep=" " secondInputIndices} )
		for i in ${!array_vcfs2[@]}; do
		ln -s ${array_indices2[i]} $(dirname ${array_vcfs2[i]})
		done

		gatk CrosscheckFingerprints -I ~{sep=" -I " firstInputs} -SI ~{sep=" -SI " secondInputs} -H ~{haplotypeDatabase} -O ~{basename}.crosscheck
	>>>

	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " HDD"
		memory: "16 GB"
	}

	output {
		File crosscheck = "~{basename}.crosscheck"
	}
}