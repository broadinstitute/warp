version 1.0

import "../../tasks/broad/RNAWithUMIsTasks.wdl" as tasks

## Copyright Broad Institute, 2021
##
## This WDL pipeline implements UMI Aware Duplicate Marking
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

workflow UMIAwareDuplicateMarking {

  input {
    File aligned_bam 
    File unaligned_bam
    String output_basename
    Boolean remove_duplicates
  }

  parameter_meta {
    aligned_bam: "Unsorted aligned bam (the output of STAR in multithread mode is not query-name sorted)"
    unaligned_bam: "Query-name sorted unaligned bam; contains UMIs in the RX tag"
    output_basename: "Basename for file outputs from this workflow"
    remove_duplicates: "If true, remove (rather than mark) duplicate reads from the output"
  }

  call tasks.SortSamByQueryName as SortSamByQueryNameAfterAlignment {
    input:
      input_bam = aligned_bam,
      output_bam_basename = output_basename + "_queryname_sorted"
  }

  call tasks.TransferReadTags {
    input:
      aligned_bam = SortSamByQueryNameAfterAlignment.output_bam,
      ubam = unaligned_bam,
      output_basename = output_basename + "_queryname_sorted_with_RX"
  }

  # First sort the aligned bam by coordinate, so we can group duplicate sets using UMIs in the next step.
  call tasks.SortSamByCoordinate as SortSamByCoordinateFirstPass {
    input:
      input_bam = TransferReadTags.output_bam,
      output_bam_basename = output_basename + ".STAR_aligned.coordinate_sorted"
  }

  # Further divide each duplicate set (a set of reads with the same insert start and end coordinates)
  # into subsets that share the same UMIs i.e. differenciate PCR duplicates from biological duplicates.
  # (biological duplicates are independent DNA molecules that are sheared such that the inserts are indistinguishable.)
  # input: a coordinate sorted bam
  # output: a coordinate sorted bam with UMIs (what are the generated tags?) .
  
  call tasks.GroupByUMIs {
    input:
      bam = SortSamByCoordinateFirstPass.output_bam,
      bam_index = SortSamByCoordinateFirstPass.output_bam_index,
      output_bam_basename = output_basename + ".grouped_by_UMI"
  }

  call tasks.SortSamByQueryName as SortSamByQueryNameBeforeDuplicateMarking {
    input:
      input_bam = GroupByUMIs.grouped_bam,
      output_bam_basename = output_basename + ".grouped.queryname_sorted"
  }

  call tasks.MarkDuplicatesUMIAware as MarkDuplicates {
    input:
      bam = SortSamByQueryNameBeforeDuplicateMarking.output_bam,
      output_basename = output_basename,
      remove_duplicates = remove_duplicates
  }

  call tasks.SortSamByCoordinate as SortSamByCoordinateSecondPass {
    input:
      input_bam = MarkDuplicates.duplicate_marked_bam,
      output_bam_basename = output_basename + ".duplicate_marked.coordinate_sorted"
  }

  output {
    File duplicate_marked_bam = SortSamByCoordinateSecondPass.output_bam
    File duplicate_marked_bam_index = SortSamByCoordinateSecondPass.output_bam_index
    File duplicate_metrics = MarkDuplicates.duplicate_metrics
  }
}
