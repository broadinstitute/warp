version 1.0

import "../../../../tasks/wdl/ImputationBeagleTasks.wdl" as beagleTasks
import "../../../../tasks/wdl/ImputationTasks.wdl" as tasks

# This workflow performs array imputation using Beagle. It's designed to scale
# to approximately 1000 samples and be used as a subworkflow for ImputationBeagle.wdl,
# which can handle larger sample sizes by splitting into batches and then merging results.

workflow ImputationBeagleBatch {
  # if this changes, update the batch_pipeline_version value in ImputationBeagle.wdl
  String pipeline_version = "1.0.0"

  input {
    Array[Array[File]] pre_chunked_multi_sample_vcfs # pre-chunked multi-sample VCFs organized by contig and chunk
    Array[Array[Int]] starts_with_overlaps # start positions with overlaps for each chunk, organized by contig and chunk
    Array[Array[Int]] ends_with_overlaps # end positions with overlaps for each chunk, organized by contig and chunk
    Array[Array[Int]] starts # start positions without overlaps for each chunk, organized by contig and chunk
    Array[Array[Int]] ends # end positions without overlaps for each chunk, organized by contig and chunk

    Boolean impute_with_allele_probabilities = false # set to true if multiple batches will be merged
    
    File ref_dict # for reheadering / adding contig lengths in the header of the output VCF
    Array[String] contigs_to_process # list of contigs that will be processed
    String reference_panel_path_prefix # path + file prefix to the bucket where the reference panel files are stored for all contigs
    String genetic_maps_path # path to the bucket where genetic maps are stored for all contigs
    String output_basename # the basename for intermediate and output files

    String? pipeline_header_line # optional additional header lines to add to the output VCF
    Float? min_dr2_for_inclusion # minimum dr2 to include a variant in the output vcf, applied after reannotation

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

  scatter (contig_index in range(length(contigs_to_process))) {
    String contig = contigs_to_process[contig_index]
    # define contig-specific filenames for the reference panel and genetic map
    String reference_basename = reference_panel_path_prefix + "." + contig
    String genetic_map_filename = genetic_maps_path + "plink." + contig + ".GRCh38.withchr.map"
    String bref3_filename = reference_basename + bref3_suffix
    String unique_variant_ids_filename = reference_basename + unique_variant_ids_suffix

    scatter (chunk_index in range(length(pre_chunked_multi_sample_vcfs[contig_index]))) {
      String chunk_basename = "${contig}_chunk_${chunk_index}"

      call beagleTasks.Phase {
        input:
          dataset_vcf = pre_chunked_multi_sample_vcfs[contig_index][chunk_index],
          ref_panel_bref3 = bref3_filename,
          chrom = contig,
          basename = "${chunk_basename}.phased",
          genetic_map_file = genetic_map_filename,
          start = starts_with_overlaps[contig_index][chunk_index],
          end = ends_with_overlaps[contig_index][chunk_index],
          cpu = beagle_cpu,
          memory_mb = beagle_phase_memory_in_gb * 1024
        }

      call beagleTasks.Impute {
        input:
        dataset_vcf = Phase.vcf,
        ref_panel_bref3 = bref3_filename,
        chrom = contig,
        basename = "${chunk_basename}.imputed",
        genetic_map_file = genetic_map_filename,
        start = starts_with_overlaps[contig_index][chunk_index],
        end = ends_with_overlaps[contig_index][chunk_index],
        impute_with_allele_probabilities = impute_with_allele_probabilities,
        cpu = beagle_cpu,
        memory_mb = beagle_impute_memory_in_gb * 1024
      }
    
      call beagleTasks.LocalizeAndSubsetVcfToRegion {
          input:
          vcf = Impute.vcf,
          start = starts[contig_index][chunk_index],
          end = ends[contig_index][chunk_index],
          contig = contig,
          output_basename = "${chunk_basename}.imputed.no_overlaps",
          gatk_docker = gatk_docker
      }

      # need to update header before gathering
      call tasks.UpdateHeader {
      input:
        vcf = LocalizeAndSubsetVcfToRegion.output_vcf,
        vcf_index = LocalizeAndSubsetVcfToRegion.output_vcf_index,
        ref_dict = ref_dict,
        basename = "${chunk_basename}.imputed.no_overlaps.update_header",
        disable_sequence_dictionary_validation = false,
        pipeline_header_line = pipeline_header_line,
        gatk_docker = gatk_docker
      }
    }

    # gather contig-wide VCFs
    call beagleTasks.GatherVcfsNoIndex as GatherVcfsNoIndexContig {
    input:
      input_vcfs = UpdateHeader.output_vcf,
      output_vcf_basename = "${output_basename}.${contig}.imputed",
      gatk_docker = gatk_docker
    }

    call beagleTasks.CreateVcfIndex as CreateIndexForGatheredVcfContig {
      input:
        vcf_input = GatherVcfsNoIndexContig.output_vcf,
        gatk_docker = gatk_docker,
        preemptible = 0
    }

  }

  output {
    Array[File] imputed_multi_sample_vcfs = CreateIndexForGatheredVcfContig.output_vcf
    Array[File] imputed_multi_sample_vcf_indexes = CreateIndexForGatheredVcfContig.output_vcf_index
  }
}
