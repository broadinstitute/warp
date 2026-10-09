version 1.0

import "../../../../tasks/wdl/ImputationTasks.wdl" as tasks
import "../../../../tasks/wdl/ImputationBeagleTasks.wdl" as beagleTasks

workflow ImputationBeagleCheckChunks {

  # if this changes, update the check_chunks_version value in ImputationBeagle.wdl
  String pipeline_version = "1.0.0"

  input {
    Int chunk_length
    Int chunk_overlaps # the padding that will be added to the beginning and end of each chunk to reduce edge effects

    File multi_sample_vcf
    File multi_sample_vcf_index

    File ref_dict # for calculating contig lengths
    Array[String] contigs_to_process # list of contigs that will be processed, based on input data
    String reference_panel_path_prefix # path + file prefix to the bucket where the reference panel files are stored for all contigs
    String output_basename # the basename for intermediate and output files

    # file extensions used to find reference panel files
    String unique_variant_ids_suffix = ".unique_variants"

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.0.0"
    String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"

    Int? error_count_override
  }

  Float chunk_length_float = chunk_length

  scatter (contig in contigs_to_process) {
    # these are specific to hg38 - contig is format 'chr1'
    String reference_basename = "${reference_panel_path_prefix}.${contig}"

    call tasks.CalculateChromosomeLength {
      input:
        ref_dict = ref_dict,
        chrom = contig,
        ubuntu_docker = ubuntu_docker
    }

    call beagleTasks.ExtractUniqueVariantIds as ExtractUniqueVariantIdsRawChromosome {
      input:
        vcf = multi_sample_vcf,
        vcf_index = multi_sample_vcf_index,
        chrom = contig
    }

    Int num_chunks = ceil(CalculateChromosomeLength.chrom_length / chunk_length_float)

    scatter (i in range(num_chunks)) {
      String chunk_contig = contig

      Int start = (i * chunk_length) + 1
      Int start_with_overlaps = if (start - chunk_overlaps < 1) then 1 else start - chunk_overlaps
      Int end = if (CalculateChromosomeLength.chrom_length < ((i + 1) * chunk_length)) then CalculateChromosomeLength.chrom_length else ((i + 1) * chunk_length)
      Int end_with_overlaps = if (CalculateChromosomeLength.chrom_length < end + chunk_overlaps) then CalculateChromosomeLength.chrom_length else end + chunk_overlaps
      String qc_scatter_position_chunk_basename = "${contig}_chunk_${i}"

      # generate the chunked vcf file that will be used for imputation, including overlaps
      call tasks.GenerateChunk {
        input:
          vcf = multi_sample_vcf,
          vcf_index = multi_sample_vcf_index,
          start = start_with_overlaps,
          end = end_with_overlaps,
          chrom = contig,
          basename = qc_scatter_position_chunk_basename,
          gatk_docker = gatk_docker
      }

      # count variants in chunk (not including overlaps) and check overlap with ref panel
      call beagleTasks.ExtractUniqueVariantIds as ExtractUniqueVariantsFilteredChunk {
        input:
        vcf = GenerateChunk.output_vcf,
        vcf_index = GenerateChunk.output_vcf_index,
        chrom = contig,
        start = start,
        end = end
      }

      call beagleTasks.CountUniqueVariantIdsInOverlap {
        input:
          variant_ids_1 = ExtractUniqueVariantsFilteredChunk.unique_variant_ids,
          variant_ids_2 = reference_basename + unique_variant_ids_suffix
      }

      call beagleTasks.CheckChunks {
        input:
          var_in_original = ExtractUniqueVariantsFilteredChunk.unique_variant_count,
          var_also_in_reference = CountUniqueVariantIdsInOverlap.var_overlap
      }
    }

    Array[File] chunked_vcfs_with_overlaps_for_imputation_by_contig = GenerateChunk.output_vcf
    Array[Int] starts_with_overlaps_by_contig = start_with_overlaps
    Array[Int] ends_with_overlaps_by_contig = end_with_overlaps
    Array[Int] starts_by_contig = start
    Array[Int] ends_by_contig = end

    call beagleTasks.CountValidContigChunks {
      input:
        valids = CheckChunks.valid
      }

    # if any chunk for any chromosome fail CheckChunks, then we fail this workflow (do not proceed to imputation)
    Int n_failed_chunks_int = select_first([error_count_override, CountValidContigChunks.n_invalid_chunks])
    call beagleTasks.ErrorWithMessageIfErrorCountNotZero as FailQCNChunks {
      input:
        errorCount = n_failed_chunks_int,
        message = "contig ${contig} had ${n_failed_chunks_int} failing chunks"
    }
  }


  call beagleTasks.StoreMetricsInfo {
    input:
      chunk_chroms = flatten(chunk_contig),
      starts = flatten(start),
      ends = flatten(end),
      vars_in_array = flatten(ExtractUniqueVariantsFilteredChunk.unique_variant_count),
      vars_in_panel = flatten(CountUniqueVariantIdsInOverlap.var_overlap),
      valids = flatten(CheckChunks.valid),
      chroms = contigs_to_process,
      vars_in_raw_input = ExtractUniqueVariantIdsRawChromosome.unique_variant_count,
      basename = output_basename
  }
  

  output {
    Array[Array[File]] chunked_vcfs_with_overlaps_for_imputation = chunked_vcfs_with_overlaps_for_imputation_by_contig
    Array[Array[Int]] starts_with_overlaps = starts_with_overlaps_by_contig
    Array[Array[Int]] ends_with_overlaps = ends_with_overlaps_by_contig
    Array[Array[Int]] starts = starts_by_contig
    Array[Array[Int]] ends = ends_by_contig
    File chunks_info = StoreMetricsInfo.chunks_info
    File failed_chunks = StoreMetricsInfo.failed_chunks
    File contigs_info = StoreMetricsInfo.contigs_info
    File n_failed_chunks = StoreMetricsInfo.n_failed_chunks
  }
}
