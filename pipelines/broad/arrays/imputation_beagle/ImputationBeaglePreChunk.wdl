version 1.0

import "../../../../tasks/broad/ImputationTasks.wdl" as tasks
import "../../../../tasks/broad/Utilities.wdl" as utils

workflow ImputationBeagle {

    String pipeline_version = "0.0.1"

    input {
        Int chunkLength = 25000000
        Int chunkOverlaps = 5000000 # this is the padding that will be added to the beginning and end of each chunk to reduce edge effects

        File multi_sample_vcf

        File ref_dict # for reheadering / adding contig lengths in the header of the ouptut VCF, and calculating contig lengths
        Array[String] contigs
        String reference_panel_path # path to the bucket where the reference panel files are stored for all contigs
        String genetic_maps_path # path to the bucket where genetic maps are stored for all contigs
        String output_basename # the basename for intermediate and output files

        Int chunks_fail_threshold = 1 # require fewer than this many chunks to fail in order to pass

        # file extensions used to find reference panel files
        String interval_list_suffix = ".interval_list"
        String bref3_suffix = ".bref3"
    }

    call tasks.CountSamples {
        input:
            vcf = multi_sample_vcf,
    }

    call tasks.PreSplitVcf {
        input:
            contigs = contigs,
            vcf = multi_sample_vcf
    }

    scatter (contig_index in range(length(contigs))) {
        # these are specific to hg38 - contig is format 'chr1'
        String reference_filename = reference_panel_path + "hgdp.tgp.gwaspy.merged." + contigs[contig_index] + ".merged.AN_added.bcf.ac2"
        String genetic_map_filename = genetic_maps_path + "plink." + contigs[contig_index] + ".GRCh38.withchr.map"

        ReferencePanelContig referencePanelContig = {
                                                        "interval_list": reference_filename + interval_list_suffix,
                                                        "bref3": reference_filename + bref3_suffix,
                                                        "contig": contigs[contig_index],
                                                        "genetic_map": genetic_map_filename
                                                    }

        call tasks.CalculateChromosomeLength {
            input:
                ref_dict = ref_dict,
                chrom = referencePanelContig.contig
        }

        call tasks.PreChunkVcf {
            input:
                chromosome_length=CalculateChromosomeLength.chrom_length,
                chunk_length = chunkLength,
                chunk_overlap = chunkOverlaps,
                chrom = contigs[contig_index],
                vcf = PreSplitVcf.chr_split_vcfs[contig_index],
                vcf_index = PreSplitVcf.chr_split_vcf_indices[contig_index]
        }

        scatter (i in range(length(PreChunkVcf.generate_chunk_vcfs))) {
            String chunk_contig = referencePanelContig.contig
            String chunk_basename = referencePanelContig.contig + "_chunk_" + i

            Int start = PreChunkVcf.starts[i]
            Int end = PreChunkVcf.ends[i]

            call tasks.CountVariantsInChunksBeagle {
                input:
                    vcf = PreChunkVcf.generate_chunk_vcfs[i],
                    vcf_index =PreChunkVcf.generate_chunk_vcf_indices[i],
                    panel_interval_list = referencePanelContig.interval_list
            }

            call tasks.CheckChunksBeagle {
                input:
                    var_in_original = CountVariantsInChunksBeagle.var_in_original,
                    var_also_in_reference = CountVariantsInChunksBeagle.var_also_in_reference
            }

            call tasks.SetIDs as SetIdsVcfToImpute {
                input:
                    vcf = PreChunkVcf.subset_vcfs[i],
                    output_basename = "input_samples_with_variant_ids"
            }

            call tasks.ExtractIDs as ExtractIdsVcfToImpute {
                input:
                    vcf = SetIdsVcfToImpute.output_vcf,
                    output_basename = "imputed_sites"
            }

            if (!CheckChunksBeagle.valid) {
                call utils.ErrorWithMessage as FailedChunk {
                    input:
                        message = chunk_basename + " was not a valid chunk"
                }
            }
            call tasks.PhaseAndImputeBeagle {
                input:
                    dataset_vcf = PreChunkVcf.generate_chunk_vcfs[i],
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
                    basename = output_basename
            }
        }

        Array[File] chromosome_vcfs = select_all(InterleaveVariants.output_vcf)
    }

    call tasks.GatherVcfs {
        input:
            input_vcfs = flatten(chromosome_vcfs),
            output_vcf_basename = output_basename + ".imputed"
    }

    call tasks.StoreChunksInfo {
        input:
            chroms = flatten(chunk_contig),
            starts = flatten(start),
            ends = flatten(end),
            vars_in_array = flatten(CountVariantsInChunksBeagle.var_in_original),
            vars_in_panel = flatten(CountVariantsInChunksBeagle.var_also_in_reference),
            valids = flatten(CheckChunksBeagle.valid),
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

struct ReferencePanelContig {
    File interval_list
    File bref3
    String contig
    File genetic_map
}
