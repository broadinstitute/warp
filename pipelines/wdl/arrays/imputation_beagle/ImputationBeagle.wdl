version 1.0

import "./ImputationBeagleCheckChunks.wdl" as ImputationBeagleCheckChunks
import "./ImputationBeagleBatch.wdl" as ImputationBeagleBatch
import "../../../../tasks/wdl/ImputationTasks.wdl" as tasks
import "../../../../tasks/wdl/ImputationBeagleTasks.wdl" as beagleTasks

workflow ImputationBeagle {
  String pipeline_version = "4.1.1"
  String check_chunks_version = "0.0.1"
  String batch_pipeline_version = "0.0.1"
  String input_qc_version = "1.3.2"
  String quota_consumed_version = "1.1.1"

  input {
    Int chunk_length = 25000000
    Int chunk_overlaps = 2000000 # the padding that will be added to the beginning and end of each chunk to reduce edge effects
    Int sample_chunk_size = 1000 # the number of samples that will be processed in parallel in each chunked scatter

    File multi_sample_vcf

    File ref_dict # for reheadering / adding contig lengths in the header of the output VCF, and calculating contig lengths
    Array[String] contigs # list of possible contigs that will be processed. note the workflow will not error out if any of these contigs are missing
    String reference_panel_path_prefix # path + file prefix to the bucket where the reference panel files are stored for all contigs
    String genetic_maps_path # path to the bucket where genetic maps are stored for all contigs
    String output_basename # the basename for intermediate and output files

    String pipeline_header_line = "" # optional additional header lines to add to the output VCF. empty string will not be added.
    Float min_dr2_for_inclusion = 0.0 # minimum dr2 to include a variant in the output vcf, applied after reannotation

    # file extensions used to find reference panel files
    String bref3_suffix = ".bref3"
    String unique_variant_ids_suffix = ".unique_variants"

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.0.0"
    String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"

    Int? error_count_override
    # the following are used to define the resources for Beagle tasks
    Int beagle_cpu = 8
    Int beagle_phase_memory_in_gb = 40
    Int beagle_impute_memory_in_gb = 45
  }

  call beagleTasks.CreateVcfIndex {
    input:
      vcf_input = multi_sample_vcf,
      gatk_docker = gatk_docker
  }

  call tasks.CountSamples {
    input:
      vcf = CreateVcfIndex.output_vcf
  }

  Float sample_chunk_size_float = sample_chunk_size
  Int num_sample_batches = ceil(CountSamples.nSamples / sample_chunk_size_float)

  call beagleTasks.CalculateContigsToProcess {
    input:
      vcf_input = CreateVcfIndex.output_vcf,
      allowed_contigs = contigs,
      gatk_docker = gatk_docker
  }

  Array[String] contigs_to_process = CalculateContigsToProcess.contigs_to_process

  call ImputationBeagleCheckChunks.ImputationBeagleCheckChunks as CheckChunks {
    input:
      chunk_length = chunk_length,
      chunk_overlaps = chunk_overlaps,
      multi_sample_vcf = CreateVcfIndex.output_vcf,
      multi_sample_vcf_index = CreateVcfIndex.output_vcf_index,
      ref_dict = ref_dict,
      contigs_to_process = contigs_to_process,
      reference_panel_path_prefix = reference_panel_path_prefix,
      genetic_maps_path = genetic_maps_path,
      output_basename = output_basename,
      unique_variant_ids_suffix = unique_variant_ids_suffix,
      gatk_docker = gatk_docker,
      ubuntu_docker = ubuntu_docker,
      error_count_override = error_count_override
  }

  # GATHER CHECK CHUNKS OUTPUT: top level is contig, next level is chunks, so [[chr1.chunk0, chr1.chunk1], [chr2.chunk0, chr2.chunk1], ...]
  Array[Array[File]] chunked_vcfs_with_overlaps_for_imputation = CheckChunks.chunked_vcfs_with_overlaps_for_imputation

  Boolean multiple_sample_batches = num_sample_batches > 1

  scatter (sample_batch_index in range(num_sample_batches)) {
    # sample FORMAT fields in vcfs start after the 8 mandatory fields plus FORMAT (FORMAT isnt mandatory
    # but if you have samples in your vcf they are).  `cut` is 1 indexed, so we start at the 10th column.
    Int batch_start_sample = (sample_batch_index * sample_chunk_size) + 10
    Int batch_end_sample = if (CountSamples.nSamples <= ((sample_batch_index + 1) * sample_chunk_size)) then CountSamples.nSamples + 9 else ((sample_batch_index + 1) * sample_chunk_size) + 9

    # only cut sample batches if there is more than one
    if (multiple_sample_batches) {
      scatter (contig_index in range(length(contigs_to_process))) {
        scatter (chunk_index in range(length(chunked_vcfs_with_overlaps_for_imputation[contig_index]))) {
          call beagleTasks.SelectSamplesWithCut {
            input:
              vcf = chunked_vcfs_with_overlaps_for_imputation[contig_index][chunk_index],
              cut_start_field = batch_start_sample,
              cut_end_field = batch_end_sample,
              basename = "filtered_input." + contigs_to_process[contig_index] + "_chunk_" + chunk_index + ".sample_batch_" + sample_batch_index
          }
        }
        Array[File] filtered_input_vcfs_for_contig = SelectSamplesWithCut.output_vcf
      }
    }

    Array[Array[File]] filtered_input_vcfs = select_first([filtered_input_vcfs_for_contig, chunked_vcfs_with_overlaps_for_imputation])

    call ImputationBeagleBatch.ImputationBeagleBatch as RunBatch {
      input:
        pre_chunked_multi_sample_vcfs = filtered_input_vcfs,
        starts_with_overlaps = CheckChunks.starts_with_overlaps,
        ends_with_overlaps = CheckChunks.ends_with_overlaps,
        starts = CheckChunks.starts,
        ends = CheckChunks.ends,
        impute_with_allele_probabilities = multiple_sample_batches,
        ref_dict = ref_dict,
        contigs_to_process = contigs_to_process,
        reference_panel_path_prefix = reference_panel_path_prefix,
        genetic_maps_path = genetic_maps_path,
        output_basename = output_basename,
        bref3_suffix = bref3_suffix,
        unique_variant_ids_suffix = unique_variant_ids_suffix,
        gatk_docker = gatk_docker,
        ubuntu_docker = ubuntu_docker,
        error_count_override = error_count_override,
        pipeline_header_line = pipeline_header_line,
        min_dr2_for_inclusion = min_dr2_for_inclusion
    }
  }

  # GATHER IMPUTATION BATCH OUTPUT [batch][chr] = File
  Array[Array[File]] imputed_multi_sample_vcf_batches = RunBatch.imputed_multi_sample_vcfs
  Array[Array[File]] imputed_multi_sample_vcf_index_batches = RunBatch.imputed_multi_sample_vcf_indexes

  # TRANSPOSE to [contig][batch] = File
  Array[Array[File]] imputed_multi_sample_vcfs_by_contig = transpose(imputed_multi_sample_vcf_batches)
  Array[Array[File]] imputed_multi_sample_vcf_indexes_by_contig = transpose(imputed_multi_sample_vcf_index_batches)

  scatter (contig_index in range(length(contigs_to_process))) {
    String contig_basename = output_basename + "." + contigs_to_process[contig_index]

    # only merge sample chunks if there is more than one
    if (multiple_sample_batches) {
      scatter (batch_index in range(length(imputed_multi_sample_vcfs_by_contig[contig_index]))) {
        call beagleTasks.QuerySampleChunkedVcfForReannotation {
          input:
            vcf = imputed_multi_sample_vcfs_by_contig[contig_index][batch_index],
          }

        call beagleTasks.RemoveAPAnnotations {
          input:
            vcf = imputed_multi_sample_vcfs_by_contig[contig_index][batch_index],
            vcf_index = imputed_multi_sample_vcf_indexes_by_contig[contig_index][batch_index],
        }

        call beagleTasks.AggregateDSandAPValuesChunked {
          input:
            query_file = QuerySampleChunkedVcfForReannotation.output_query_file,
            n_samples = QuerySampleChunkedVcfForReannotation.n_samples,
        }
      }

      call beagleTasks.MergeSampleChunksVcfsWithPaste {
        input:
          input_vcfs = RemoveAPAnnotations.output_vcf,
          output_vcf_basename = contig_basename + ".imputed.samples_merged",
      }

      call beagleTasks.CreateVcfIndex as IndexMergedSampleChunksVcfs {
        input:
          vcf_input = MergeSampleChunksVcfsWithPaste.output_vcf,
          gatk_docker = gatk_docker
      }
    
      call beagleTasks.AggregateChunkedDR2AndAF {
        input:
          sample_chunked_annotation_files = AggregateDSandAPValuesChunked.output_summary_file
      }

      call beagleTasks.ReannotateDR2AndAF {
        input:
        vcf = IndexMergedSampleChunksVcfs.output_vcf,
        vcf_index = IndexMergedSampleChunksVcfs.output_vcf_index,
        annotations_tsv = AggregateChunkedDR2AndAF.output_annotations_file,
        annotations_tsv_index = AggregateChunkedDR2AndAF.output_annotations_file_index
      }
    }

    # Define contig VCF for all input samples
    File all_samples_contig_vcf = select_first([ReannotateDR2AndAF.output_vcf, imputed_multi_sample_vcf_batches[0][contig_index]])
    File all_samples_contig_vcf_index = select_first([ReannotateDR2AndAF.output_vcf_index, imputed_multi_sample_vcf_index_batches[0][contig_index]])

    # only filter by dr2 if the user has defined a threshold greater than 0
    if (min_dr2_for_inclusion > 0.0) {
      call beagleTasks.FilterVcfByDR2 {
          input:
          vcf = all_samples_contig_vcf,
          vcf_index = all_samples_contig_vcf_index,
          basename = contig_basename + ".imputed",
          dr2_threshold = min_dr2_for_inclusion,
          gatk_docker = gatk_docker
      }
    }

    # Define contig VCF for all input samples with filtering
    File all_samples_contig_vcf_filtered = select_first([FilterVcfByDR2.output_vcf, all_samples_contig_vcf])
    File all_samples_contig_vcf_index_filtered = select_first([FilterVcfByDR2.output_vcf_index, all_samples_contig_vcf_index])
  }
  
  output {
    Array[File] imputed_multi_sample_vcfs = all_samples_contig_vcf_filtered
    Array[File] imputed_multi_sample_vcf_indexes = all_samples_contig_vcf_index_filtered
    File chunks_info = CheckChunks.chunks_info
    File contigs_info = CheckChunks.contigs_info
  }

  meta {
    allowNestedInputs: true
  }
}
