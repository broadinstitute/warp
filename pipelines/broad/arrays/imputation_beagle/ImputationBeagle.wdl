version 1.0

import "../../../../structs/imputation/ImputationBeagleStructs.wdl" as structs
import "../../../../tasks/broad/ImputationTasks.wdl" as tasks
import "../../../../tasks/broad/ImputationBeagleTasks.wdl" as beagleTasks

workflow ImputationBeagle {

  String pipeline_version = "1.0.0"

  input {
    Int chunkLength = 25000000
    Int chunkOverlaps = 2000000 # this is the padding that will be added to the beginning and end of each chunk to reduce edge effects

    File multi_sample_vcf

    File ref_dict # for reheadering / adding contig lengths in the header of the ouptut VCF, and calculating contig lengths
    Array[String] contigs
    String reference_panel_path_prefix # path + file prefix to the bucket where the reference panel files are stored for all contigs
    String genetic_maps_path # path to the bucket where genetic maps are stored for all contigs
    String output_basename # the basename for intermediate and output files

    # file extensions used to find reference panel files
    String bed_suffix = ".bed"
    String bref3_suffix = ".bref3"

    String gatk_docker = "broadinstitute/gatk:4.6.0.0"
    String ubuntu_docker = "ubuntu:20.04"

    Int? error_count_override
  }

  call tasks.CountSamples {
    input:
      vcf = multi_sample_vcf
  }

  call tasks.CreateVcfIndex {
    input:
      vcf_input = multi_sample_vcf,
      gatk_docker = gatk_docker
  }

  Float chunkLengthFloat = chunkLength

  scatter (contig in contigs) {
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
      String chunk_basename = referencePanelContig.contig + "_chunk_" + i

      # generate the chunked vcf file that will be used for imputation, including overlaps
      call tasks.GenerateChunk {
        input:
          vcf = CreateVcfIndex.vcf,
          vcf_index = CreateVcfIndex.vcf_index,
          start = startWithOverlaps,
          end = endWithOverlaps,
          chrom = referencePanelContig.contig,
          basename = chunk_basename,
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

    scatter (i in range(num_chunks)) {
      String chunk_basename_imputed = referencePanelContig.contig + "_chunk_" + i + "_imputed"

      # max amount of cpus you can ask for is 96 so at a max of 10k samples we can only ask for 9 cpu a sample.
      # these values are based on trying to optimize for pre-emptibility using a 400k sample referene panel
      # and up to a 10k sample input vcf
      Int beagle_cpu = if (CountSamples.nSamples <= 1000) then 8 else floor(CountSamples.nSamples / 1000) * 9
      Int beagle_phase_memory_in_gb = if (CountSamples.nSamples <= 1000) then 22 else ceil(beagle_cpu * 1.5)
      Int beagle_impute_memory_in_gb = if (CountSamples.nSamples <= 1000) then 30 else ceil(beagle_cpu * 4.3)

      # max amount of cpus you can ask for is 96 so at a max of 10k samples we can only ask for 9 cpu a sample.
      # these values are based on trying to optimize for pre-emptibility using a 400k sample referene panel
      # and up to a 10k sample input vcf
      Int beagle_cpu = if (CountSamples.nSamples <= 1000) then 8 else floor(CountSamples.nSamples / 1000) * 9
      Int beagle_phase_memory_in_gb = if (CountSamples.nSamples <= 1000) then 22 else ceil(beagle_cpu * 1.5)
      Int beagle_impute_memory_in_gb = if (CountSamples.nSamples <= 1000) then 30 else ceil(beagle_cpu * 4.3)

      call beagleTasks.Phase {
        input:
          dataset_vcf = chunkedVcfsWithOverlapsForImputation[i],
          ref_panel_bref3 = referencePanelContig.bref3,
          chrom = referencePanelContig.contig,
          basename = chunk_basename_imputed,
          genetic_map_file = referencePanelContig.genetic_map,
          start = start[i],
          end = end[i],
          cpu = beagle_cpu,
          memory_mb = beagle_phase_memory_in_gb * 1024
      }

      call beagleTasks.Impute {
        input:
          dataset_vcf = Phase.vcf,
          ref_panel_bref3 = referencePanelContig.bref3,
          chrom = referencePanelContig.contig,
          basename = chunk_basename_imputed,
          genetic_map_file = referencePanelContig.genetic_map,
          start = start[i],
          end = end[i],
          cpu = beagle_cpu,
          memory_mb = beagle_impute_memory_in_gb * 1024
      }

      call tasks.CreateVcfIndex as IndexImputedBeagle {
        input:
          vcf_input = Impute.vcf,
          gatk_docker = gatk_docker
      }

      call tasks.UpdateHeader {
        input:
          vcf = IndexImputedBeagle.vcf,
          vcf_index = IndexImputedBeagle.vcf_index,
          ref_dict = ref_dict,
          basename = chunk_basename_imputed,
          disable_sequence_dictionary_validation = false,
          gatk_docker = gatk_docker
      }

      call tasks.SeparateMultiallelics {
        input:
          original_vcf = UpdateHeader.output_vcf,
          original_vcf_index = UpdateHeader.output_vcf_index,
          output_basename = chunk_basename_imputed
      }

      call tasks.RemoveSymbolicAlleles {
        input:
          original_vcf = SeparateMultiallelics.output_vcf,
          original_vcf_index = SeparateMultiallelics.output_vcf_index,
          output_basename = chunk_basename_imputed,
          gatk_docker = gatk_docker
      }
    }

    Array[File] chromosome_vcfs = select_all(RemoveSymbolicAlleles.output_vcf)
  }

  call tasks.GatherVcfs {
    input:
      input_vcfs = flatten(chromosome_vcfs),
      output_vcf_basename = output_basename + ".imputed",
      gatk_docker = gatk_docker
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
    File imputed_multi_sample_vcf = GatherVcfs.output_vcf
    File imputed_multi_sample_vcf_index = GatherVcfs.output_vcf_index
    File chunks_info = StoreChunksInfo.chunks_info
  }

  meta {
    allowNestedInputs: true
  }

}
