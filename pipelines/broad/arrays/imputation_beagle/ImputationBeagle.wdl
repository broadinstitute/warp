version 1.0

import "../../../../tasks/broad/ImputationTasks.wdl" as tasks
import "../../../../tasks/broad/Utilities.wdl" as utils

workflow ImputationBeagle {

  String pipeline_version = "0.0.1"

  input {
    Int chunkLength = 25000000
    Int chunkOverlaps = 5000000 # this is the padding that will be added to the beginning and end of each chunk to reduce edge effects

    File multi_sample_vcf
    File multi_sample_vcf_index

    Boolean perform_extra_qc_steps = false # these are optional additional extra QC steps from Amit's group that should only be
    # run for large sample sets, especially a diverse set of samples (it's further limiting called at sites to 95% and by HWE)
    Float? optional_qc_max_missing
    Float? optional_qc_hwe
    File ref_dict # for reheadering / adding contig lengths in the header of the ouptut VCF, and calculating contig lengths
    Array[String] contigs
    String reference_panel_path # path to the bucket where the reference panel files are stored for all contigs
    String genetic_maps_path # path to the bucket where genetic maps are stored for all contigs
    String output_callset_name # the output callset name
    Boolean split_output_to_single_sample = false
    Int merge_ssvcf_mem_mb = 3000 # the memory allocation for MergeSingleSampleVcfs (in mb)

    Float frac_above_maf_5_percent_well_imputed_threshold = 0.9 # require fraction of maf > 0.05 sites well imputed to be greater than this to pass
    Int chunks_fail_threshold = 1 # require fewer than this many chunks to fail in order to pass

    # file extensions used to find reference panel files
    String vcf_suffix = ".vcf.gz"
    String vcf_index_suffix = ".vcf.gz.tbi"
    String bref3_suffix = ".bref3"
  }

  call tasks.CountSamples {
    input:
      vcf = multi_sample_vcf,
  }

  Float chunkLengthFloat = chunkLength

  scatter (contig in contigs) {
    # these are specific to hg38 - contig is format 'chr1'
    String reference_filename = reference_panel_path + "hgdp.tgp.gwaspy.merged." + contig + ".merged.AN_added.bcf.ac2"
    String genetic_map_filename = genetic_maps_path + "plink." + contig + ".GRCh38.withchr.map"

    ReferencePanelContig referencePanelContig = {
      "vcf": reference_filename + vcf_suffix,
      "vcf_index": reference_filename + vcf_index_suffix,
      "bref3": reference_filename + bref3_suffix,
      "contig": contig,
      "genetic_map": genetic_map_filename
    }

    call tasks.CalculateChromosomeLength {
      input:
        ref_dict = ref_dict,
        chrom = referencePanelContig.contig
    }

    call tasks.CreateRefPanelIntervalLists {
      input:
        ref_panel_vcf = referencePanelContig.vcf,
        ref_panel_vcf_index = referencePanelContig.vcf_index
    }

    Int num_chunks = ceil(CalculateChromosomeLength.chrom_length / chunkLengthFloat)

    scatter (i in range(num_chunks)) {
      String chunk_contig = referencePanelContig.contig
      Int start = (i * chunkLength) + 1
      Int startWithOverlaps = if (start - chunkOverlaps < 1) then 1 else start - chunkOverlaps
      Int end = if (CalculateChromosomeLength.chrom_length < ((i + 1) * chunkLength)) then CalculateChromosomeLength.chrom_length else ((i + 1) * chunkLength)
      Int endWithOverlaps = if (CalculateChromosomeLength.chrom_length < end + chunkOverlaps) then CalculateChromosomeLength.chrom_length else end + chunkOverlaps
      String chunk_basename = referencePanelContig.contig + "_chunk_" + i

      call tasks.GenerateChunk {
        input:
          vcf = multi_sample_vcf,
          vcf_index = multi_sample_vcf_index,
          start = startWithOverlaps,
          end = endWithOverlaps,
          chrom = referencePanelContig.contig,
          basename = chunk_basename
      }

      if (perform_extra_qc_steps) {
        call tasks.OptionalQCSites {
          input:
            input_vcf = GenerateChunk.output_vcf,
            input_vcf_index = GenerateChunk.output_vcf_index,
            output_vcf_basename = chunk_basename,
            optional_qc_max_missing = optional_qc_max_missing,
            optional_qc_hwe = optional_qc_hwe
        }
      }

      call tasks.CountVariantsInChunksBeagle {
        input:
          vcf = select_first([OptionalQCSites.output_vcf,  GenerateChunk.output_vcf]),
          vcf_index = select_first([OptionalQCSites.output_vcf_index, GenerateChunk.output_vcf_index]),
          panel_interval_list = CreateRefPanelIntervalLists.interval_list
      }

      call tasks.CheckChunksBeagle {
        input:
          var_in_original = CountVariantsInChunksBeagle.var_in_original,
          var_in_reference = CountVariantsInChunksBeagle.var_in_reference
      }

      call tasks.SubsetVcfToRegion {
        input:
          vcf = multi_sample_vcf,
          vcf_index = multi_sample_vcf_index,
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

      if (CheckChunksBeagle.valid) {
        call tasks.PhaseAndImputeBeagle {
          input:
            dataset_vcf = select_first([OptionalQCSites.output_vcf,  GenerateChunk.output_vcf]),
            ref_panel_bref3 = referencePanelContig.bref3,
            chrom = referencePanelContig.contig,
            basename = chunk_basename,
            genetic_map_file = referencePanelContig.genetic_map,
            start = start,
            end = end
        }

        call tasks.UpdateHeader {
          input:
            vcf = PhaseAndImputeBeagle.vcf,
            vcf_index = PhaseAndImputeBeagle.vcf_index,
            ref_dict = ref_dict,
            basename = chunk_basename + "_imputed"
        }

        call tasks.SeparateMultiallelics {
          input:
            original_vcf = UpdateHeader.output_vcf,
            original_vcf_index = UpdateHeader.output_vcf_index,
            output_basename = chunk_basename + "_imputed"
        }

        call tasks.RemoveSymbolicAlleles {
          input:
            original_vcf = SeparateMultiallelics.output_vcf,
            original_vcf_index = SeparateMultiallelics.output_vcf_index,
            output_basename = chunk_basename + "_imputed"
        }

        call tasks.SetIDs {
          input:
            vcf = RemoveSymbolicAlleles.output_vcf,
            output_basename = chunk_basename + "_imputed"
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

    Array[File] chromosome_vcfs = select_all(InterleaveVariants.output_vcf)
  }

  call tasks.GatherVcfs {
    input:
      input_vcfs = flatten(chromosome_vcfs),
      output_vcf_basename = output_callset_name + ".imputed"
  }

  call tasks.StoreChunksInfo {
    input:
      chroms = flatten(chunk_contig),
      starts = flatten(start),
      ends = flatten(end),
      vars_in_array = flatten(CountVariantsInChunksBeagle.var_in_original),
      vars_in_panel = flatten(CountVariantsInChunksBeagle.var_in_reference),
      valids = flatten(CheckChunksBeagle.valid),
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
    File imputed_multi_sample_vcf = GatherVcfs.output_vcf
    File imputed_multi_sample_vcf_index = GatherVcfs.output_vcf_index
    File chunks_info = StoreChunksInfo.chunks_info
    File failed_chunks = StoreChunksInfo.failed_chunks
    File n_failed_chunks = StoreChunksInfo.n_failed_chunks
  }

  meta {
    allowNestedInputs: true
  }

}

struct ReferencePanelContig {
  File vcf
  File vcf_index
  File bref3
  String contig
  File genetic_map
}
