version 1.0

import "../../../../structs/imputation/ImputationStructs.wdl" as structs
import "../../../../tasks/wdl/ImputationTasks.wdl" as tasks
import "../../../../tasks/wdl/Utilities.wdl" as utils

workflow Imputation {

  String pipeline_version = "1.1.22"

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

    Boolean perform_extra_qc_steps = false # these are optional additional extra QC steps from Amit's group that should only be
    # run for large sample sets, especially a diverse set of samples (it's further limiting called at sites to 95% and by HWE)
    Float? optional_qc_max_missing
    Float? optional_qc_hwe
    File ref_dict # for reheadering / adding contig lengths in the header of the ouptut VCF, and calculating contig lengths
    Array[String] contigs
    String reference_panel_path # path to the bucket where the reference panel files are stored for all contigs
    File genetic_maps_eagle
    String output_callset_name # the output callset name
    Boolean split_output_to_single_sample = false
    Int merge_ssvcf_mem_mb = 3000 # the memory allocation for MergeSingleSampleVcfs (in mb)

    Float frac_above_maf_5_percent_well_imputed_threshold = 0.9 # require fraction of maf > 0.05 sites well imputed to be greater than this to pass
    Int chunks_fail_threshold = 1 # require fewer than this many chunks to fail in order to pass

    # file extensions used to find reference panel files
    String vcf_suffix = ".vcf.gz"
    String vcf_index_suffix = ".vcf.gz.tbi"
    String bcf_suffix = ".bcf"
    String bcf_index_suffix =  ".bcf.csi"
    String m3vcf_suffix = ".cleaned.m3vcf.gz"
  }

  if (defined(single_sample_vcfs) && defined(multi_sample_vcf)) {
    call utils.ErrorWithMessage as ErrorMessageDoubleInput{
      input:
        message = "single_sample_vcfs and multi_sample_vcf cannot both be defined as input"
    }
  }

  if (!defined(single_sample_vcfs) && !defined(multi_sample_vcf)) {
    call utils.ErrorWithMessage as ErrorMessageNoInput {
      input:
        message = "One (and only one) of single_sample_vcfs and multi_sample_vcf must be defined as input"
    }
  }

  if (defined(single_sample_vcfs)) {
    call tasks.MergeSingleSampleVcfs {
      input:
        input_vcfs = select_first([single_sample_vcfs]),
        input_vcf_indices = select_first([single_sample_vcf_indices]),
        output_vcf_basename = "merged_input_samples",
        memory_mb = merge_ssvcf_mem_mb
    }
  }

  File vcf_to_impute = select_first([multi_sample_vcf, MergeSingleSampleVcfs.output_vcf])
  File vcf_index_to_impute = select_first([multi_sample_vcf_index, MergeSingleSampleVcfs.output_vcf_index])

  call tasks.CountSamples {
    input:
      vcf = vcf_to_impute,
  }


  Float chunkLengthFloat = chunkLength

  scatter (contig in contigs) {

    String reference_filename = reference_panel_path + "ALL.chr" + contig + ".phase3_integrated.20130502.genotypes.cleaned"

    ReferencePanelContig referencePanelContig = {
      "vcf": reference_filename + vcf_suffix,
      "vcf_index": reference_filename + vcf_index_suffix,
      "bcf": reference_filename + bcf_suffix,
      "bcf_index": reference_filename + bcf_index_suffix,
      "m3vcf": reference_filename + m3vcf_suffix,
      "contig": contig
    }

    call tasks.CalculateChromosomeLength {
      input:
        ref_dict = ref_dict,
        chrom = referencePanelContig.contig
    }

    Int num_chunks = ceil(CalculateChromosomeLength.chrom_length / chunkLengthFloat)

    scatter (i in range(num_chunks)) {
      String chunk_contig = referencePanelContig.contig
      Int start = (i * chunkLength) + 1
      Int startWithOverlaps = if (start - chunkOverlaps < 1) then 1 else start - chunkOverlaps
      Int end = if (CalculateChromosomeLength.chrom_length < ((i + 1) * chunkLength)) then CalculateChromosomeLength.chrom_length else ((i + 1) * chunkLength)
      Int endWithOverlaps = if (CalculateChromosomeLength.chrom_length < end + chunkOverlaps) then CalculateChromosomeLength.chrom_length else end + chunkOverlaps

      call tasks.GenerateChunk {
        input:
          vcf = vcf_to_impute,
          vcf_index = vcf_index_to_impute,
          start = startWithOverlaps,
          end = endWithOverlaps,
          chrom = referencePanelContig.contig,
          basename = "chrom_" + referencePanelContig.contig + "_chunk_" + i
      }

      if (perform_extra_qc_steps) {
        call tasks.OptionalQCSites {
          input:
            input_vcf = GenerateChunk.output_vcf,
            input_vcf_index = GenerateChunk.output_vcf_index,
            output_vcf_basename =  "chrom_" + referencePanelContig.contig + "_chunk_" + i,
            optional_qc_max_missing = optional_qc_max_missing,
            optional_qc_hwe = optional_qc_hwe
        }
      }

      call tasks.CountVariantsInChunks {
        input:
          vcf = select_first([OptionalQCSites.output_vcf,  GenerateChunk.output_vcf]),
          vcf_index = select_first([OptionalQCSites.output_vcf_index, GenerateChunk.output_vcf_index]),
          panel_vcf = referencePanelContig.vcf,
          panel_vcf_index = referencePanelContig.vcf_index
      }
      call tasks.CheckChunks {
        input:
          vcf = select_first([OptionalQCSites.output_vcf,  GenerateChunk.output_vcf]),
          vcf_index = select_first([OptionalQCSites.output_vcf_index, GenerateChunk.output_vcf_index]),
          panel_vcf = referencePanelContig.vcf,
          panel_vcf_index = referencePanelContig.vcf_index,
          var_in_original = CountVariantsInChunks.var_in_original,
          var_in_reference = CountVariantsInChunks.var_in_reference
      }

      call tasks.SubsetVcfToRegion {
        input:
          vcf = vcf_to_impute,
          vcf_index = vcf_index_to_impute,
          output_basename = "input_samples_subset_to_chunk",
          contig = referencePanelContig.contig,
          start = start,
          end = end
      }

      call tasks.SetIDs as SetIdsVcfToImpute {
        input:
          vcf = SubsetVcfToRegion.output_vcf,
          output_basename = "input_samples_with_variant_ids"
      }

      call tasks.ExtractIDs as ExtractIdsVcfToImpute {
        input:
          vcf = SetIdsVcfToImpute.output_vcf,
          output_basename = "imputed_sites"
      }

      if (CheckChunks.valid) {
        call tasks.PhaseVariantsEagle {
          input:
            dataset_bcf = CheckChunks.valid_chunk_bcf,
            dataset_bcf_index = CheckChunks.valid_chunk_bcf_index,
            reference_panel_bcf = referencePanelContig.bcf,
            reference_panel_bcf_index = referencePanelContig.bcf_index,
            chrom = referencePanelContig.contig,
            genetic_map_file = genetic_maps_eagle,
            start = startWithOverlaps,
            end = endWithOverlaps
        }

        call tasks.Minimac4 {
          input:
            ref_panel = referencePanelContig.m3vcf,
            phased_vcf = PhaseVariantsEagle.dataset_prephased_vcf,
            prefix = "chrom_" + referencePanelContig.contig + "_chunk_" + i +"_imputed",
            chrom = referencePanelContig.contig,
            start = start,
            end = end,
            window = chunkOverlaps
        }

        call tasks.AggregateImputationQCMetrics {
          input:
            infoFile = Minimac4.info,
            nSamples = CountSamples.nSamples,
            basename = output_callset_name + "chrom_" + referencePanelContig.contig + "_chunk_" + i
        }

        call tasks.UpdateHeader {
          input:
            vcf = Minimac4.vcf,
            vcf_index = Minimac4.vcf_index,
            ref_dict = ref_dict,
            basename = "chrom_" + referencePanelContig.contig + "_chunk_" + i +"_imputed"
        }

        call tasks.SeparateMultiallelics {
          input:
            original_vcf = UpdateHeader.output_vcf,
            original_vcf_index = UpdateHeader.output_vcf_index,
            output_basename = "chrom" + referencePanelContig.contig + "_chunk_" + i +"_imputed"
        }

        call tasks.RemoveSymbolicAlleles {
          input:
            original_vcf = SeparateMultiallelics.output_vcf,
            original_vcf_index = SeparateMultiallelics.output_vcf_index,
            output_basename = "chrom" + referencePanelContig.contig + "_chunk_" + i +"_imputed"
        }

        call tasks.SetIDs {
          input:
            vcf = RemoveSymbolicAlleles.output_vcf,
            output_basename = "chrom" + referencePanelContig.contig + "_chunk_" + i +"_imputed"
        }

        call tasks.ExtractIDs {
          input:
            vcf = SetIDs.output_vcf,
            output_basename = "imputed_sites"
        }
      }
      call tasks.FindSitesUniqueToFileTwoOnly {
        input:
          file1 = select_first([ExtractIDs.ids, write_lines([])]),
          file2 = ExtractIdsVcfToImpute.ids
      }

      call tasks.SelectVariantsByIds {
        input:
          vcf = SetIdsVcfToImpute.output_vcf,
          vcf_index = SetIdsVcfToImpute.output_vcf_index,
          ids = FindSitesUniqueToFileTwoOnly.missing_sites,
          basename = "imputed_sites_to_recover"
      }

      call tasks.RemoveAnnotations {
        input:
          vcf = SelectVariantsByIds.output_vcf,
          basename = "imputed_sites_to_recover_annotations_removed"
      }

      call tasks.InterleaveVariants {
        input:
          vcfs = select_all([RemoveAnnotations.output_vcf, SetIDs.output_vcf]),
          basename = output_callset_name
      }
    }
    Array[File] aggregatedImputationMetrics = select_all(AggregateImputationQCMetrics.aggregated_metrics)
    Array[File] chromosome_vcfs = select_all(InterleaveVariants.output_vcf)
  }

  Array[String] phased_vcfs = flatten(chromosome_vcfs)

  call tasks.GetMissingContigList {
    input:
      ref_dict = ref_dict,
      included_contigs = write_lines(contigs)
  }

  scatter (missing_contig in GetMissingContigList.missing_contigs) {
    call tasks.CalculateChromosomeLength as CalculateMissingChromosomeLength {
      input:
        ref_dict = ref_dict,
        chrom = missing_contig
    }

    Int num_chunks_missing_contig = ceil(CalculateMissingChromosomeLength.chrom_length / chunkLengthFloat)

    scatter (i_missing_contig in range(num_chunks_missing_contig)) {
      Int start_missing_contig = (i_missing_contig * chunkLength) + 1
      Int end_missing_contig = if (CalculateMissingChromosomeLength.chrom_length < ((i_missing_contig + 1) * chunkLength)) then CalculateMissingChromosomeLength.chrom_length else ((i_missing_contig + 1) * chunkLength)

      call tasks.SubsetVcfToRegion as SubsetVcfToRegionMissingContig{
        input:
          vcf = vcf_to_impute,
          vcf_index = vcf_index_to_impute,
          output_basename = "input_samples_subset_to_chunk",
          contig = missing_contig,
          start = start_missing_contig,
          end = end_missing_contig,
          exclude_filtered = true
      }

      call tasks.SetIDs as SetIDsMissingContigs {
        input:
          vcf = SubsetVcfToRegionMissingContig.output_vcf,
          output_basename = "unimputed_contigs_" + missing_contig +"_"+ i_missing_contig + "_with_ids"
      }

      call tasks.RemoveAnnotations as RemoveAnnotationsMissingContigs {
        input:
          vcf = SetIDsMissingContigs.output_vcf,
          basename = "unimputed_contigs_" + missing_contig +"_"+ i_missing_contig + "_annotations_removed"
      }
    }
  }

  Array[String] missing_remove_annotation_vcfs = flatten(RemoveAnnotationsMissingContigs.output_vcf)

  scatter(missing_remove_annotation_vcf in missing_remove_annotation_vcfs){
    call tasks.ReplaceHeader {
      input:
        vcf_to_replace_header = missing_remove_annotation_vcf,
        vcf_with_new_header = phased_vcfs[0]
    }
  }

  Array[String] missing_contig_vcfs = ReplaceHeader.output_vcf
  Array[String] unsorted_vcfs = flatten([phased_vcfs, missing_contig_vcfs])

  call tasks.GatherVcfs {
    input:
      input_vcfs = unsorted_vcfs,
      output_vcf_basename = output_callset_name
  }

  call tasks.MergeImputationQCMetrics {
    input:
      metrics = flatten(aggregatedImputationMetrics),
      basename = output_callset_name
  }

  if (MergeImputationQCMetrics.frac_above_maf_5_percent_well_imputed < frac_above_maf_5_percent_well_imputed_threshold) {
    call utils.ErrorWithMessage as FailQCWellImputedFrac {
      input:
        message = "Well imputed fraction was " + MergeImputationQCMetrics.frac_above_maf_5_percent_well_imputed + ", QC failure threshold was set at " + frac_above_maf_5_percent_well_imputed_threshold
    }
  }

  call tasks.StoreChunksInfo {
    input:
      chroms = flatten(chunk_contig),
      starts = flatten(start),
      ends = flatten(end),
      vars_in_array = flatten(CountVariantsInChunks.var_in_original),
      vars_in_panel = flatten(CountVariantsInChunks.var_in_reference),
      valids = flatten(CheckChunks.valid),
      basename = output_callset_name
  }
  
  Int n_failed_chunks_int = read_int(StoreChunksInfo.n_failed_chunks)

  if (n_failed_chunks_int >= chunks_fail_threshold) {
    call utils.ErrorWithMessage as FailQCNChunks {
      input:
        message = n_failed_chunks_int + " chunks failed imputation, QC threshold was set to " + chunks_fail_threshold
    }
  }

  if (split_output_to_single_sample) {
    call tasks.SplitMultiSampleVcf {
      input:
        multiSampleVcf = GatherVcfs.output_vcf,
        nSamples = CountSamples.nSamples
    }
  }


  output {
    Array[File]? imputed_single_sample_vcfs = SplitMultiSampleVcf.single_sample_vcfs
    Array[File]? imputed_single_sample_vcf_indices = SplitMultiSampleVcf.single_sample_vcf_indices
    File imputed_multisample_vcf = GatherVcfs.output_vcf
    File imputed_multisample_vcf_index = GatherVcfs.output_vcf_index
    File aggregated_imputation_metrics = MergeImputationQCMetrics.aggregated_metrics
    File chunks_info = StoreChunksInfo.chunks_info
    File failed_chunks = StoreChunksInfo.failed_chunks
    File n_failed_chunks = StoreChunksInfo.n_failed_chunks
  }

  meta {
    allowNestedInputs: true
  }

}
