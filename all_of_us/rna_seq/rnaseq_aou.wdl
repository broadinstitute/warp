# import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:samtofastq_v1-0_BETA/versions/8/plain-WDL/descriptor" as samtofastq_wdl
# import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:star_v1-0_BETA/versions/8/plain-WDL/descriptor" as star_wdl
# import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:markduplicates_v1-0_BETA/versions/6/plain-WDL/descriptor" as markduplicates_wdl
# import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:rsem_v1-0_BETA/versions/6/plain-WDL/descriptor" as rsem_wdl
# import "https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_gtex:rnaseqc2_v1-0_BETA/versions/4/plain-WDL/descriptor" as rnaseqc_wdl

import "./samtofastq.wdl" as samtofastq_wdl
import "./star.wdl" as star_wdl
import "./markduplicates.wdl" as markduplicates_wdl
import "./rsem.wdl" as rsem_wdl
import "./rnaseqc2.wdl" as rnaseqc_wdl


workflow rnaseq_pipeline_bam_workflow {

    String prefix

    call samtofastq_wdl.samtofastq {
        input: prefix=prefix
    }

    call star_align_wdl.star {
        input: fastq1=samtofastq.fastq1, fastq2=samtofastq.fastq2, prefix=prefix
    }

    call prersem_wdl.remove_IDS_reads {
        input: transcriptome_bam=star.transcriptome_bam, prefix=prefix
    }

    call rsem_wdl.rsem {
        input: transcriptome_bam=remove_IDS_reads.transcriptome_noIDS_bam, prefix=prefix
    }

#    call bamsync_wdl.bamsync {
#        input: target_bam=star.bam_file, target_bam_index=star.bam_index, prefix=prefix
#    }

    call markduplicates_wdl.markduplicates {
        #input: input_bam=bamsync.patched_bam_file, prefix=prefix
        input: input_bam=star.bam_file, prefix=prefix

    }

    call rnaseqc_wdl.rnaseqc2 {
        input: bam_file=markduplicates.bam_file, sample_id=prefix
    }
}