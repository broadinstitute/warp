version 1.0

# excised from https://github.com/broadinstitute/warp/blob/848b3a7c2d9b9e1b525e8276809f29d2b60e63af/tasks/broad/Utilities.wdl
workflow BamToCram {

  input {
    File    input_bam
    File    ref_fasta
    File    ref_fasta_index
    String  output_basename
    Int     preemptible_tries   =   3
  }

  call ConvertToCram {
    input:
      input_bam         = input_bam,
      ref_fasta         = ref_fasta,
      ref_fasta_index   = ref_fasta_index,
      output_basename   = output_basename,
      preemptible_tries =   preemptible_tries
  }

  output {
     File   output_cram         =    ConvertToCram.output_cram
     File   output_cram_index   =   ConvertToCram.output_cram_index
     File   output_cram_md5     =   ConvertToCram.output_cram_md5
  }

  meta {
    allowNestedInputs: true
  }
  
}

# Convert BAM file to CRAM format
# Note that reading CRAMs directly with Picard is not yet supported
task ConvertToCram {
  input {
    File input_bam
    File ref_fasta
    File ref_fasta_index
    String output_basename
    Int preemptible_tries = 3

    Int disk_size = ceil((2 * size(input_bam, "GiB")) + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")) + 20
  }

  command <<<
    set -e
    set -o pipefail

    samtools view -C -T ~{ref_fasta} ~{input_bam} | \
    tee ~{output_basename}.cram | \
    md5sum | awk '{print $1}' > ~{output_basename}.cram.md5

    # Create REF_CACHE. Used when indexing a CRAM
    seq_cache_populate.pl -root ./ref/cache ~{ref_fasta}
    export REF_PATH=:
    export REF_CACHE=./ref/cache/%2s/%2s/%s

    samtools index ~{output_basename}.cram
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
    preemptible: preemptible_tries
    memory: "3 GiB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_cram = "~{output_basename}.cram"
    File output_cram_index = "~{output_basename}.cram.crai"
    File output_cram_md5 = "~{output_basename}.cram.md5"
  }
}