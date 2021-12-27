version 1.0

import "./RNAWithUMIsTasks.wdl" as tasks


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

  String pipeline_version = "0.1.0"

  input {
    File aligned_bam 
    String output_basename
  }

  parameter_meta {
    aligned_bam: "Aligned bam sorted by the query (read) name"
    output_basename: "Basename for file outputs from this workflow"
  }

  # First sort the aligned bam by coordinate, so we can group duplicate sets using UMIs in the next step.
  call tasks.SortSamUMIAware as SortSamFirst {
    input:
      input_bam = aligned_bam,
      output_bam_basename = output_basename + ".STAR_aligned.coorinate_sorted",
      sort_order = "coordinate"
  }

  # Further divide each duplicate set (a set of reads with the same insert start and end coordinates)
  # into subsets that share the same UMIs i.e. differenciate PCR duplicates from biological duplicates.
  # (biological duplicates are independent DNA molecules that are sheared such that the inserts are indistinguishable.)
  # input: a coordinate sorted bam
  # output: a coordinate sorted bam with UMIs (what are the generated tags?) .
  
  call tasks.GroupByUMIs {
    input:
      bam = SortSamFirst.output_bam,
      bam_index = select_first([SortSamFirst.output_bam_index, "bam_index_not_found"]),
      output_bam_basename = output_basename + ".grouped_by_UMI"
  }

  call tasks.SortSamUMIAware as SortSamQueryName {
    input:
      input_bam = GroupByUMIs.grouped_bam,
      output_bam_basename = output_basename + ".grouped.queryname_sorted",
      sort_order = "queryname"
  }

  call tasks.MarkDuplicatesUMIAware as MarkDuplicates {
    input:
      bam = SortSamQueryName.output_bam,
      output_basename = output_basename
  }

  call tasks.SortSamUMIAware as SortSamSecond {
    input:
      input_bam = MarkDuplicates.duplicate_marked_bam,
      output_bam_basename = output_basename + ".duplicate_marked.coordinate_sorted",
      sort_order = "coordinate"
  }

  output {
    File duplicate_marked_bam = SortSamSecond.output_bam
    File duplicate_marked_bam_index = select_first([SortSamSecond.output_bam_index, "bam_index_not_found"])
    File duplicate_metrics = MarkDuplicates.duplicate_metrics
  }
}
