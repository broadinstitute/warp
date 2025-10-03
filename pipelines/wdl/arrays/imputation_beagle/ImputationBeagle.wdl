version 1.0

import "../../../../structs/imputation/ImputationBeagleStructs.wdl" as structs
import "../../../../tasks/wdl/ImputationTasks.wdl" as tasks
import "../../../../tasks/wdl/ImputationBeagleTasks.wdl" as beagleTasks

workflow ImputationBeagle {
  String pipeline_version = "2.2.1"
  String input_qc_version = "1.2.1"
  String quota_consumed_version = "1.1.0"

  input {
    Int chunkLength = 25000000
    Int chunkOverlaps = 2000000 # this is the padding that will be added to the beginning and end of each chunk to reduce edge effects
    Int sample_chunk_size = 1000 # this is the number of samples that will be processed in parallel in each chunked scatter

    File multi_sample_vcf

    File ref_dict # for reheadering / adding contig lengths in the header of the ouptut VCF, and calculating contig lengths
    Array[String] contigs # list of possible contigs that will be processed. note the workflow will not error out if any of these contigs are missing
    String reference_panel_path_prefix # path + file prefix to the bucket where the reference panel files are stored for all contigs
    String genetic_maps_path # path to the bucket where genetic maps are stored for all contigs
    String output_basename # the basename for intermediate and output files

    String? pipeline_header_line # optional additional header lines to add to the output VCF
    Float? min_dr2_for_inclusion # minimum dr2 to include a variant in the output vcf, applied after reannotation

    # file extensions used to find reference panel files
    String bed_suffix = ".bed"
    String bref3_suffix = ".bref3"

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.0.0"
    String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"

    Int? error_count_override
    # the following are used to define the resources for Beagle tasks
    Int beagle_cpu = 8
    Int beagle_phase_memory_in_gb = 40
    Int beagle_impute_memory_in_gb = 45
  }

  # have to define these here to use them in nested scatters
  String defined_pipeline_header_line = if defined(pipeline_header_line) then select_first([pipeline_header_line]) else ""
  Float defined_min_dr2_for_inclusion = if defined(min_dr2_for_inclusion) then select_first([min_dr2_for_inclusion]) else 0.0

  call tasks.CountSamples {
    input:
      vcf = multi_sample_vcf
  }

  Float sample_chunk_size_float = sample_chunk_size
  Int num_sample_chunks = ceil(CountSamples.nSamples / sample_chunk_size_float)

  call beagleTasks.CreateVcfIndex {
    input:
      vcf_input = multi_sample_vcf,
      gatk_docker = gatk_docker
  }

  call beagleTasks.CalculateContigsToProcess {
    input:
      vcf_input = multi_sample_vcf,
      allowed_contigs = contigs,
      gatk_docker = gatk_docker
  }

  Array[String] contigs_to_process = CalculateContigsToProcess.contigs_to_process

  Float chunkLengthFloat = chunkLength

  scatter (contig in contigs_to_process) {
    # these are specific to hg38 - contig is format 'chr1'
    String reference_basename = reference_panel_path_prefix + "." + contig
    String genetic_map_filename = genetic_maps_path + "plink." + contig + ".GRCh38.withchr.map"

    ReferencePanelContig referencePanelContig = {
      "bed": reference_basename + bed_suffix,
      "bref3": reference_basename  + bref3_suffix,
      "contig": contig,
      "genetic_map": genetic_map_filename
    }


    call tasks.CalculateChromosomeLength {
      input:
        ref_dict = ref_dict,
        chrom = referencePanelContig.contig,
        ubuntu_docker = ubuntu_docker
    }

    Int num_chunks = ceil(CalculateChromosomeLength.chrom_length / chunkLengthFloat)

    scatter (i in range(num_chunks)) {
      String chunk_contig = referencePanelContig.contig

      Int start = (i * chunkLength) + 1
      Int startWithOverlaps = if (start - chunkOverlaps < 1) then 1 else start - chunkOverlaps
      Int end = if (CalculateChromosomeLength.chrom_length < ((i + 1) * chunkLength)) then CalculateChromosomeLength.chrom_length else ((i + 1) * chunkLength)
      Int endWithOverlaps = if (CalculateChromosomeLength.chrom_length < end + chunkOverlaps) then CalculateChromosomeLength.chrom_length else end + chunkOverlaps
      String qc_scatter_position_chunk_basename = referencePanelContig.contig + "_chunk_" + i

      # generate the chunked vcf file that will be used for imputation, including overlaps
      call tasks.GenerateChunk {
        input:
          vcf = CreateVcfIndex.output_vcf,
          vcf_index = CreateVcfIndex.output_vcf_index,
          start = startWithOverlaps,
          end = endWithOverlaps,
          chrom = referencePanelContig.contig,
          basename = qc_scatter_position_chunk_basename,
          gatk_docker = gatk_docker
      }

      call beagleTasks.CountVariantsInChunks {
        input:
          vcf = GenerateChunk.output_vcf,
          vcf_index = GenerateChunk.output_vcf_index,
          panel_bed_file = referencePanelContig.bed,
          gatk_docker = gatk_docker
      }

      call beagleTasks.CheckChunks {
        input:
          var_in_original = CountVariantsInChunks.var_in_original,
          var_also_in_reference = CountVariantsInChunks.var_also_in_reference
      }
    }

    Array[File] chunkedVcfsWithOverlapsForImputation = GenerateChunk.output_vcf

    call tasks.StoreChunksInfo as StoreContigLevelChunksInfo {
      input:
        chroms = chunk_contig,
        starts = start,
        ends = end,
        vars_in_array = CountVariantsInChunks.var_in_original,
        vars_in_panel = CountVariantsInChunks.var_also_in_reference,
        valids = CheckChunks.valid,
        basename = output_basename
    }

    # if any chunk for any chromosome fail CheckChunks, then we will not impute run any task in the next scatter,
    # namely phasing and imputing which would be the most costly to throw away
    Int n_failed_chunks_int = select_first([error_count_override, read_int(StoreContigLevelChunksInfo.n_failed_chunks)])
    call beagleTasks.ErrorWithMessageIfErrorCountNotZero as FailQCNChunks {
      input:
        errorCount = n_failed_chunks_int,
        message = "contig " + referencePanelContig.contig + " had " + n_failed_chunks_int + " failing chunks"
    }
  }

  scatter (contig_index in range(length(contigs_to_process))) {
    # cant have the same variable names in different scatters, so add _2 to "differentiate"
    String reference_basename_2 = reference_panel_path_prefix + "." + contigs_to_process[contig_index]
    String genetic_map_filename_2 = genetic_maps_path + "plink." + contigs_to_process[contig_index] + ".GRCh38.withchr.map"

    ReferencePanelContig referencePanelContig_2 = {
                                                  "bed": reference_basename_2 + bed_suffix,
                                                  "bref3": reference_basename_2  + bref3_suffix,
                                                  "contig": contigs_to_process[contig_index],
                                                  "genetic_map": genetic_map_filename_2
                                                }
    Int num_chunks_2 = num_chunks[contig_index]
    scatter (i in range(num_chunks_2)) {
      String impute_scatter_position_chunk_basename = referencePanelContig_2.contig + "_chunk_" + i

      scatter (j in range(num_sample_chunks)) {
        # sample FORMAT fields in vcfs start after the 8 mandatory fields plus FORMAT (FORMAT isnt mandatory
        # but if you have samples in your vcf they are).  `cut` is 1 indexed, so we start at the 10th column.
        Int start_sample = (j * sample_chunk_size) + 10
        Int end_sample = if (CountSamples.nSamples <= ((j + 1) * sample_chunk_size)) then CountSamples.nSamples + 9 else ((j + 1) * sample_chunk_size ) + 9
        String impute_scatter_sample_chunk_basename = impute_scatter_position_chunk_basename + ".sample_chunk_" + j
        Boolean impute_with_allele_probablities = num_sample_chunks > 1

        # only cut sample chunks if there is more than one
        if (num_sample_chunks > 1) {
          call beagleTasks.SelectSamplesWithCut {
            input:
              vcf = chunkedVcfsWithOverlapsForImputation[contig_index][i],
              cut_start_field = start_sample,
              cut_end_field = end_sample,
              basename = impute_scatter_sample_chunk_basename
          }
        }

        call beagleTasks.Phase {
          input:
            dataset_vcf = select_first([SelectSamplesWithCut.output_vcf, chunkedVcfsWithOverlapsForImputation[contig_index][i]]),
            ref_panel_bref3 = referencePanelContig_2.bref3,
            chrom = referencePanelContig_2.contig,
            basename = impute_scatter_sample_chunk_basename + ".phased",
            genetic_map_file = referencePanelContig_2.genetic_map,
            start = startWithOverlaps[contig_index][i],
            end = endWithOverlaps[contig_index][i],
            cpu = beagle_cpu,
            memory_mb = beagle_phase_memory_in_gb * 1024,
            for_dependency = FailQCNChunks.done
        }

        call beagleTasks.Impute {
          input:
            dataset_vcf = Phase.vcf,
            ref_panel_bref3 = referencePanelContig_2.bref3,
            chrom = referencePanelContig_2.contig,
            basename = impute_scatter_sample_chunk_basename + ".imputed",
            genetic_map_file = referencePanelContig_2.genetic_map,
            start = startWithOverlaps[contig_index][i],
            end = endWithOverlaps[contig_index][i],
            impute_with_allele_probabilities = impute_with_allele_probablities,
            cpu = beagle_cpu,
            memory_mb = beagle_impute_memory_in_gb * 1024
        }

        call beagleTasks.LocalizeAndSubsetVcfToRegion {
          input:
            vcf = Impute.vcf,
            start = start[contig_index][i],
            end = end[contig_index][i],
            contig = referencePanelContig_2.contig,
            output_basename = impute_scatter_sample_chunk_basename + ".imputed.no_overlaps",
            gatk_docker = gatk_docker
        }
        # set up DR2 and AF reannotation for the imputed VCFs
        if (num_sample_chunks > 1) {
          call beagleTasks.QuerySampleChunkedVcfForReannotation {
            input:
              vcf = LocalizeAndSubsetVcfToRegion.output_vcf,
          }

          call beagleTasks.RemoveAPAnnotations {
            input:
              vcf = LocalizeAndSubsetVcfToRegion.output_vcf,
              vcf_index = LocalizeAndSubsetVcfToRegion.output_vcf_index,
          }

          call beagleTasks.RecalculateDR2AndAFChunked {
            input:
              query_file = QuerySampleChunkedVcfForReannotation.output_query_file,
              n_samples = QuerySampleChunkedVcfForReannotation.n_samples,
          }
        }
        # create a non optional File for use in future tasks, wdl quirk requires the dummy value to exist to validate
        File ap_annotations_removed_vcf = select_first([RemoveAPAnnotations.output_vcf, "gs://fake/will_fail.txt"])
        File chunked_dr2_af = select_first([RecalculateDR2AndAFChunked.output_summary_file, "gs://fake/will_fail.txt"])
      }

      # only merge sample chunks if there is more than one
      if (num_sample_chunks > 1) {
        call beagleTasks.MergeSampleChunksVcfsWithPaste {
          input:
            input_vcfs = ap_annotations_removed_vcf,
            output_vcf_basename = impute_scatter_position_chunk_basename + ".imputed.no_overlaps.samples_merged",
        }

        call beagleTasks.CreateVcfIndex as IndexMergedSampleChunksVcfs {
          input:
            vcf_input = MergeSampleChunksVcfsWithPaste.output_vcf,
            gatk_docker = gatk_docker
        }

        call beagleTasks.AggregateChunkedDR2AndAF {
          input:
            sample_chunked_annotation_files = chunked_dr2_af
        }

        call beagleTasks.ReannotateDR2AndAF {
          input:
            vcf = IndexMergedSampleChunksVcfs.output_vcf,
            vcf_index = IndexMergedSampleChunksVcfs.output_vcf_index,
            annotations_tsv = AggregateChunkedDR2AndAF.output_annotations_file,
            annotations_tsv_index = AggregateChunkedDR2AndAF.output_annotations_file_index
        }
      }

      # only filter by dr2 if the user has defined a threshold greater than 0
      if (defined_min_dr2_for_inclusion > 0.0) {
        call beagleTasks.FilterVcfByDR2 {
          input:
            vcf = select_first([ReannotateDR2AndAF.output_vcf, LocalizeAndSubsetVcfToRegion.output_vcf[0]]),
            vcf_index = select_first([ReannotateDR2AndAF.output_vcf_index, LocalizeAndSubsetVcfToRegion.output_vcf_index[0]]),
            basename = impute_scatter_position_chunk_basename + ".imputed.no_overlaps.filtered",
            dr2_threshold = defined_min_dr2_for_inclusion,
            gatk_docker = gatk_docker
        }
      }

      call tasks.UpdateHeader {
        input:
          vcf = select_first([FilterVcfByDR2.output_vcf, ReannotateDR2AndAF.output_vcf, LocalizeAndSubsetVcfToRegion.output_vcf[0]]),
          vcf_index = select_first([FilterVcfByDR2.output_vcf_index, ReannotateDR2AndAF.output_vcf_index, LocalizeAndSubsetVcfToRegion.output_vcf_index[0]]),
          ref_dict = ref_dict,
          basename = impute_scatter_position_chunk_basename + ".imputed.no_overlaps.update_header",
          disable_sequence_dictionary_validation = false,
          pipeline_header_line = defined_pipeline_header_line,
          gatk_docker = gatk_docker
      }
    }

    Array[File] chromosome_vcfs = select_all(UpdateHeader.output_vcf)
  }

  call beagleTasks.GatherVcfsNoIndex {
    input:
      input_vcfs = flatten(chromosome_vcfs),
      output_vcf_basename = output_basename + ".imputed",
      gatk_docker = gatk_docker
  }

  call beagleTasks.CreateVcfIndex as CreateIndexForGatheredVcf {
    input:
      vcf_input = GatherVcfsNoIndex.output_vcf,
      gatk_docker = gatk_docker,
      preemptible = 0
  }

  call tasks.StoreChunksInfo {
    input:
      chroms = flatten(chunk_contig),
      starts = flatten(start),
      ends = flatten(end),
      vars_in_array = flatten(CountVariantsInChunks.var_in_original),
      vars_in_panel = flatten(CountVariantsInChunks.var_also_in_reference),
      valids = flatten(CheckChunks.valid),
      basename = output_basename
  }
  
  output {
    File imputed_multi_sample_vcf = CreateIndexForGatheredVcf.output_vcf
    File imputed_multi_sample_vcf_index = CreateIndexForGatheredVcf.output_vcf_index
    File chunks_info = StoreChunksInfo.chunks_info
  }

  meta {
    allowNestedInputs: true
  }

}
