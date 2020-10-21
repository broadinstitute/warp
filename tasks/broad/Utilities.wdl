version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines utility tasks used for processing of sequencing data.
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

# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  input {
    File ref_dict
    Int preemptible_tries
  }
  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("~{ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
    memory: "2 GiB"
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

# This task calls picard's IntervalListTools to scatter the input interval list into scatter_count sub interval lists
# Note that the number of sub interval lists may not be exactly equal to scatter_count.  There may be slightly more or less.
# Thus we have the block of python to count the number of generated sub interval lists.
task ScatterIntervalList {
  input {
    File interval_list
    Int scatter_count
    Int break_bands_at_multiples_of
  }

  command <<<
    set -e
    mkdir out
    java -Xms1g -jar /usr/gitc/picard.jar \
      IntervalListTools \
      SCATTER_COUNT=~{scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=~{break_bands_at_multiples_of} \
      INPUT=~{interval_list} \
      OUTPUT=out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals))
    CODE
  >>>
  output {
    Array[File] out = glob("out/*/*.interval_list")
    Int interval_count = read_int(stdout())
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    memory: "2 GiB"
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
    Int preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")
  Int disk_size = ceil(2 * size(input_bam, "GiB") + ref_size) + 20

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
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
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

# Convert CRAM file to BAM format
task ConvertToBam {
  input {
    File input_cram
    File ref_fasta
    File ref_fasta_index
    String output_basename
  }

  command <<<
    set -e
    set -o pipefail

    samtools view -b -o ~{output_basename}.bam -T ~{ref_fasta} ~{input_cram}

    samtools index ~{output_basename}.bam
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    preemptible: 3
    memory: "3 GiB"
    cpu: "1"
    disks: "local-disk 200 HDD"
  }
  output {
    File output_bam = "~{output_basename}.bam"
    File output_bam_index = "~{output_basename}.bam.bai"
  }
}

# Calculates sum of a list of floats
task SumFloats {
  input {
    Array[Float] sizes
    Int preemptible_tries
  }

  command <<<
    python -c "print ~{sep="+" sizes}"
  >>>
  output {
    Float total_size = read_float(stdout())
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
    preemptible: preemptible_tries
  }
}
