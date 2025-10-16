version 1.0

## Copyright Broad Institute, 2021
##
## This WDL defines tasks used for alignment of human whole-genome or exome sequencing data using Illumina's DRAGEN open source mapper.
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
task SamToFastqAndDragmapAndMba {
  input {
    File input_bam
    String output_bam_basename

    ReferenceFasta reference_fasta
    DragmapReference dragmap_reference

    Int compression_level
    Int preemptible_tries
    Boolean hard_clip_reads = false
    Boolean unmap_contaminant_reads = true

    String docker = "us.gcr.io/broad-gotc-prod/dragmap:1.1.2-1.2.1-2.26.10-1.11-1643839530"
    Int cpu = 16
    Float disk_multiplier = 8
    Int memory_mb = 40960
  }

  Float unmapped_bam_size = size(input_bam, "GiB")
  Float ref_size = size(reference_fasta.ref_fasta, "GiB") + size(reference_fasta.ref_fasta_index, "GiB") + size(reference_fasta.ref_dict, "GiB")
  Float bwa_ref_size = ref_size + size(reference_fasta.ref_alt, "GiB") + size(reference_fasta.ref_amb, "GiB") + size(reference_fasta.ref_ann, "GiB") + size(reference_fasta.ref_bwt, "GiB") + size(reference_fasta.ref_pac, "GiB") + size(reference_fasta.ref_sa, "GiB")
  Float dragmap_ref_size = size(dragmap_reference.reference_bin, "GiB") + size(dragmap_reference.hash_table_cfg_bin, "GiB") + size(dragmap_reference.hash_table_cmp, "GiB")
  Int disk_size_gb = ceil(unmapped_bam_size + bwa_ref_size + dragmap_ref_size + (disk_multiplier * unmapped_bam_size) + 20)

  command <<<
    set -euxo pipefail

    DRAGMAP_VERSION=$(dragen-os --version)

    if [ -z ${DRAGMAP_VERSION} ]; then
        exit 1;
    fi

    mkdir dragen_reference
    mv ~{dragmap_reference.reference_bin} ~{dragmap_reference.hash_table_cfg_bin} ~{dragmap_reference.hash_table_cmp} dragen_reference

    dragen-os -b ~{input_bam} -r dragen_reference --interleaved=1 --preserve-map-align-order true 2> >(tee ~{output_bam_basename}.dragmap.stderr.log >&2) | samtools view -h -O BAM - > aligned.bam
    java -Dsamjdk.compression_level=~{compression_level} -Xms1000m -Xmx1000m -jar /picard/picard.jar \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ATTRIBUTES_TO_REMOVE=RG \
      ATTRIBUTES_TO_REMOVE=NM \
      ATTRIBUTES_TO_REMOVE=MD \
      ALIGNED_BAM=aligned.bam \
      UNMAPPED_BAM=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      REFERENCE_SEQUENCE=~{reference_fasta.ref_fasta} \
      PAIRED_RUN=true \
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
      PROGRAM_RECORD_ID="dragen-os" \
      PROGRAM_GROUP_VERSION="${DRAGMAP_VERSION}" \
      PROGRAM_GROUP_COMMAND_LINE="dragen-os -b ~{input_bam} -r dragen_reference --interleaved=1" \
      PROGRAM_GROUP_NAME="dragen-os" \
      UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
      ALIGNER_PROPER_PAIR_FLAGS=true \
      UNMAP_CONTAMINANT_READS=~{unmap_contaminant_reads} \
      ADD_PG_TAG_TO_READS=false
  >>>
  runtime {
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"
    cpu: cpu
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File dragmap_stderr_log = "~{output_bam_basename}.dragmap.stderr.log"
  }
}
