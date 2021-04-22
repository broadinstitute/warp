version 1.0

import "../../../tasks/skylab/StarAlign.wdl" as StarAlignBam
import "../../../tasks/skylab/Picard.wdl" as Picard
import "../../../tasks/skylab/TrimAdapters.wdl" as TrimAdapters


#import "../../../tasks/skylab/GroupMetricsOutputs.wdl" as GroupQCs
#import "../../../tasks/skylab/LoomUtils.wdl" as LoomUtils

workflow SmartSeq2SingleNuclei {
  meta {
    description: "Process SmartSeq2 snRNA-Seq data, include reads alignment, QC metrics collection, and gene expression quantitication"
  }

  input {
    # load annotation
    #File genome_ref_fasta
    File gene_ref_flat
    # load index
    # ref index name

    File tar_star_reference
    File annotations_gtf

    # samples

    String stranded
    String input_id
    String output_name

    File adapter_list
    File fastq1
    File? fastq2
    Boolean paired_end
  }
  
  # version of this pipeline
  String pipeline_version = "1.0.0"

  parameter_meta {
    genome_ref_fasta: "Genome reference in fasta format"
    stranded: "Library strand information example values: FR RF NONE"
    input_id: "Sample name or cell_names"
    output_name: "Output name, can include path"
    fastq1: "R1 in paired end reads"
    fastq2: "R2 in paired end reads"
    paired_end: "Boolean flag denoting if the sample is paired end or not"
    tar_star_reference: "star genome reference"
    annotations_gtf: "gtf containing annotations for gene tagging (must match star reference)"
  }


     call TrimApapters.TrimAdapters {
       input:
         fastq1 = fastq1,
         fastq2 = select_first([fastq2]),
         adapter_list = adapter_list
     }



  output {
    # version of this pipeline
    String pipeline_version_out = pipeline_version

    # quality control outputs
    File aligned_bam = HISAT2_output_bam

    # data outputs
    # loom
  }
}
