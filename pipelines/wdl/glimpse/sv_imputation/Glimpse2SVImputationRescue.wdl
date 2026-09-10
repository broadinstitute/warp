version 1.0

import "./PreprocessPLsGVCF.wdl" as PreprocessPLsGVCF
import "./Glimpse2SVImputationBatch.wdl" as Glimpse2SVImputationBatch
import "../../../../tasks/wdl/Glimpse2SVImputationTasks.wdl" as Glimpse2SVImputationTasks

workflow Glimpse2SVImputation {
    String pipeline_version = "0.0.28"
    String preprocess_pls_gvcf_pipeline_version = "0.0.16"
    String batch_pipeline_version = "0.0.20"
    String quota_consumed_version = "0.0.2"
    String input_qc_version = "0.0.2"

    input {
        String output_basename

        Array[String] chromosomes

        # per chromosome inputs
        Array[File] merged_popped_contig_vcfs       # output of MergePoppedContigVcfs or RecomputePoppedAfInfo
        Array[Array[File]] annotations              # if non-empty for a contig, RecomputePoppedAfInfo will be run

        Array[Int] num_samples

        # Optional filter: variants with INFO score below this threshold will be excluded from the final output VCFs
        Float info_filter_for_inclusion = 0.0

        String merge_docker = "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    }

    scatter (contig_idx in range(length(chromosomes))) {
        File merged_popped_contig_vcf = merged_popped_contig_vcfs[contig_idx]
        Array[File] annotations_for_contig = annotations[contig_idx]

        if (length(annotations_for_contig) > 0) {
            call Glimpse2SVImputationTasks.RecomputeAndAnnotate as RecomputePoppedAfInfo {
                input:
                    merged_vcf_or_bcf = merged_popped_contig_vcf,
                    annotations = annotations_for_contig,
                    num_samples = num_samples,
                    output_basename = output_basename + "." + chromosomes[contig_idx] + ".glimpse2.popped.merged.reannotated",
                    docker_merge = merge_docker
            }
        }

        File final_popped_contig_vcf = select_first([RecomputePoppedAfInfo.merged_imputed_vcf, merged_popped_contig_vcf])

        if (info_filter_for_inclusion > 0.0) {
            call Glimpse2SVImputationTasks.FilterVcfByInfo as FilterPoppedContigByInfo {
                input:
                    vcf_or_bcf = final_popped_contig_vcf,
                    info_threshold = info_filter_for_inclusion,
                    output_basename = output_basename + "." + chromosomes[contig_idx] + ".glimpse2.popped.info_filtered"
            }
        }

        File final_filtered_popped_contig_vcf = select_first([FilterPoppedContigByInfo.output_vcf, final_popped_contig_vcf])

        call Glimpse2SVImputationTasks.CreateVcfIndexAndMd5 as IndexFinalPoppedContig {
            input:
                vcf_input_or_bcf = final_filtered_popped_contig_vcf,
                output_basename = output_basename + "." + chromosomes[contig_idx],
                gatk_docker = gatk_docker,
                preemptible = 0
        }
    }

    output {
        Array[File] imputed_vcfs = IndexFinalPoppedContig.output_vcf
        Array[File] imputed_vcf_indexes = IndexFinalPoppedContig.output_vcf_index
    }
}

