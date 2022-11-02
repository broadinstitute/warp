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
task SamToFastqAndBwaMemAndMba {
  input {
    File input_bam
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

  Float unmapped_bam_size = size(input_bam, "GiB")
  Float ref_size = size(reference_fasta.ref_fasta, "GiB") + size(reference_fasta.ref_fasta_index, "GiB") + size(reference_fasta.ref_dict, "GiB")
  Float bwa_ref_size = ref_size + size(reference_fasta.ref_alt, "GiB") + size(reference_fasta.ref_amb, "GiB") + size(reference_fasta.ref_ann, "GiB") + size(reference_fasta.ref_bwt, "GiB") + size(reference_fasta.ref_pac, "GiB") + size(reference_fasta.ref_sa, "GiB")
  # Sometimes the output is larger than the input, or a task can spill to disk.
  # In these cases we need to account for the input (1) and the output (1.5) or the input(1), the output(1), and spillage (.5).
  Float disk_multiplier = 2.5
  Int disk_size = ceil(unmapped_bam_size + bwa_ref_size + (disk_multiplier * unmapped_bam_size) + 20)

  command <<<


    # This is done before "set -o pipefail" because "bwa" will have a rc=1 and we don't want to allow rc=1 to succeed
    # because the sed may also fail with that error and that is something we actually want to fail on.
    BWA_VERSION=$(/usr/gitc/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //')

    set -o pipefail
    set -e

    if [ -z ${BWA_VERSION} ]; then
        exit 1;
    fi

    # set the bash variable needed for the command-line
    bash_ref_fasta=~{reference_fasta.ref_fasta}
    # if reference_fasta.ref_alt has data in it or allow_empty_ref_alt is set
    if [ -s ~{reference_fasta.ref_alt} ] || ~{allow_empty_ref_alt}; then
      java -Xms1000m -Xmx1000m -jar /usr/gitc/picard.jar \
        SamToFastq \
        INPUT=~{input_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      /usr/gitc/~{bwa_commandline} /dev/stdin - 2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
      java -Dsamjdk.compression_level=~{compression_level} -Xms1000m -Xmx1000m -jar /usr/gitc/picard.jar \
        MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=~{input_bam} \
        OUTPUT=~{output_bam_basename}.bam \
        REFERENCE_SEQUENCE=~{reference_fasta.ref_fasta} \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        ~{true='CLIP_OVERLAPPING_READS=true' false="" hard_clip_reads} \
        ~{true='CLIP_OVERLAPPING_READS_OPERATOR=H' false="" hard_clip_reads} \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        PROGRAM_RECORD_ID="bwamem" \
        PROGRAM_GROUP_VERSION="${BWA_VERSION}" \
        PROGRAM_GROUP_COMMAND_LINE="~{bwa_commandline}" \
        PROGRAM_GROUP_NAME="bwamem" \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=~{unmap_contaminant_reads} \
        ADD_PG_TAG_TO_READS=false

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

task SamSplitter {
  input {
    File input_bam
    Int n_reads
    Int compression_level
    Int preemptible_tries = 3
  }

  Float unmapped_bam_size = size(input_bam, "GiB")
  # Since the output bams are less compressed than the input bam we need a disk multiplier that's larger than 2.
  Float disk_multiplier = 2.5
  Int disk_size = ceil(disk_multiplier * unmapped_bam_size + 20)

  command {
    set -e
    mkdir output_dir

    total_reads=$(samtools view -c ~{input_bam})

    java -Dsamjdk.compression_level=~{compression_level} -Xms3000m -Xmx3600m -jar /usr/gitc/picard.jar SplitSamByNumberOfReads \
      INPUT=~{input_bam} \
      OUTPUT=output_dir \
      SPLIT_TO_N_READS=~{n_reads} \
      TOTAL_READS_IN_INPUT=$total_reads
  }
  output {
    Array[File] split_bams = glob("output_dir/*.bam")
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    preemptible: preemptible_tries
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
}
