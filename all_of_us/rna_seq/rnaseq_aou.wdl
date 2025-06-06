# import "https://api.firecloud.org/ga4gh/v1/tools/broad_drc_aou_aux:samToFastq/versions/1/plain-WDL/descriptor" as samtofastq_wdl
# import "https://api.firecloud.org/ga4gh/v1/tools/gtex_v10_pg:star_align_v10pg/versions/12/plain-WDL/descriptor" as star_align_wdl
# import "https://api.firecloud.org/ga4gh/v1/tools/gtex_v10_pg:pre_RSEM_processing_v10pg/versions/6/plain-WDL/descriptor" as prersem_wdl
# import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:rsem_v1-0_BETA/versions/6/plain-WDL/descriptor" as rsem_wdl
# #import "https://api.firecloud.org/ga4gh/v1/tools/broad_drc_aou_aux:bamsync_v1-0_BETA_041824_WB/versions/5/plain-WDL/descriptor" as bamsync_wdl
# import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:markduplicates_v1-0_BETA/versions/6/plain-WDL/descriptor" as markduplicates_wdl
# import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:rnaseqc2_v1-0_BETA/versions/4/plain-WDL/descriptor" as rnaseqc_wdl

import "./samtofastq.wdl" as samtofastq_wdl
import "./star.wdl" as star_align_wdl
import "./markduplicates.wdl" as markduplicates_wdl
import "./rsem.wdl" as rsem_wdl
import "./rnaseqc2.wdl" as rnaseqc_wdl


workflow rnaseq_pipeline_bam_workflow {

    String prefix
    String pipeline_version = "aou_9.0.1"

    call samtofastq_wdl.samtofastq {
        input: prefix=prefix
    }

    call star_align_wdl.star {
        input: fastq1=samtofastq.fastq1, fastq2=samtofastq.fastq2, prefix=prefix
    }

    call rsem_wdl.rsem {
        input: transcriptome_bam=star.transcriptome_bam, prefix=prefix
    }

    call markduplicates_wdl.markduplicates {
        #input: input_bam=bamsync.patched_bam_file, prefix=prefix
        input: input_bam=star.bam_file, prefix=prefix

    }

    call rnaseqc_wdl.rnaseqc2 {
        input: bam_file=markduplicates.bam_file, sample_id=prefix
    }
}