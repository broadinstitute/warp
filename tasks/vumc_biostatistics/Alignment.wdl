version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for alignment of human whole-genome or exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "../../structs/dna_seq/DNASeqStructs.wdl"

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment, then stream to MergeBamAlignment
task FastqToBwaMemAndMba {
  input {
    File fastq_1
    File fastq_2
    String bwa_commandline
    String output_bam_basename

    # reference_fasta.ref_alt is the .alt file from bwa-kit
    # (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative".
    ReferenceFasta reference_fasta

    Int compression_level
    Int preemptible_tries
    Boolean hard_clip_reads = false
    Boolean unmap_contaminant_reads = true
    Boolean allow_empty_ref_alt = false
  }

  Float fastq_size = size(fastq_1, "GiB") + size(fastq_2, "GiB")
  Float ref_size = size(reference_fasta.ref_fasta, "GiB") + size(reference_fasta.ref_fasta_index, "GiB") + size(reference_fasta.ref_dict, "GiB")
  Float bwa_ref_size = ref_size + size(reference_fasta.ref_alt, "GiB") + size(reference_fasta.ref_amb, "GiB") + size(reference_fasta.ref_ann, "GiB") + size(reference_fasta.ref_bwt, "GiB") + size(reference_fasta.ref_pac, "GiB") + size(reference_fasta.ref_sa, "GiB")
  # Sometimes the output is larger than the input, or a task can spill to disk.
  # In these cases we need to account for the input (1) and the output (1.5) or the input(1), the output(1), and spillage (.5).
  Float disk_multiplier = 2.5
  Int disk_size = ceil(fastq_size + bwa_ref_size + (disk_multiplier * fastq_size) + 20)

  command <<<

    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=~{reference_fasta.ref_fasta}
    # if reference_fasta.ref_alt has data in it or allow_empty_ref_alt is set
    if [ -s ~{reference_fasta.ref_alt} ] || ~{allow_empty_ref_alt}; then
      /usr/gitc/~{bwa_commandline} -R "@RG\tID:1\tPU:~{output_bam_basename}\tLB:~{output_bam_basename}\tSM:~{output_bam_basename}\tPL:ILLUMINA" ~{fastq_1} ~{fastq_2} - 2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
      samtools view -bS -o ~{output_bam_basename}.unsorted.bam

      if ~{!allow_empty_ref_alt}; then
        grep -m1 "read .* ALT contigs" ~{output_bam_basename}.bwa.stderr.log | \
        grep -v "read 0 ALT contigs"
      fi

    # else reference_fasta.ref_alt is empty or could not be found
    else
      echo ref_alt input is empty or not provided. >&2
      exit 1;
    fi
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    preemptible: preemptible_tries
    memory: "14 GiB"
    cpu: "16"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
  }
}

task FastqSplitter {
  input {
    File fastq
    Array[String] out_files
    Int compression_level = 1
    Int preemptible_tries = 3
  }

  Float fastq_size = size(fastq, "GiB")
  # Since the output bams are less compressed than the input bam we need a disk multiplier that's larger than 2.
  Float disk_multiplier = 2.5
  Int disk_size = ceil(disk_multiplier * fastq_size + 20)

  command {
    set -e

    fastqsplitter -i ~{fastq} -c ~{compression_level} -o ~{sep=" -o " out_files}
  }
  output {
    Array[File] split_fastqs = out_files
  }
  runtime {
    docker: "shengqh/cqs_exomeseq:latest"
    preemptible: preemptible_tries
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
}
