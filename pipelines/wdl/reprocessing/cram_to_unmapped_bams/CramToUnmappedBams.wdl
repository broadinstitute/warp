version 1.0

# Exactly one of input_cram and input_bam should be supplied to this workflow. If an input_cram is supplied, a ref_fasta
# and ref_fasta_index must also be supplied. The ref_fasta and ref_fasta_index are used to generate a bam, so if an
# input_cram is not supplied, the input_bam is used instead and the ref_fasta and ref_fasta_index are not needed.

# If the output_map file is provided, it is expected to be a tab-separated file containing a list of all the read group ids
# found in the input_cram / input_bam and the desired name of the unmapped bams generated for each.
# If the file is not provided, the output names of the unmapped bams will be the read_group_id<unmapped_bam_suffix>
workflow CramToUnmappedBams {

  String pipeline_version = "1.1.3"

  input {
    File? input_cram
    File? input_bam
    File? ref_fasta
    File? ref_fasta_index
    File? output_map
    String base_file_name
    String unmapped_bam_suffix = ".unmapped.bam"
    Int additional_disk = 20
  }

  if (defined(input_cram)) {
    Float cram_size = size(input_cram, "GiB")
    String bam_from_cram_name = basename(input_cram_path, ".cram")
    String input_cram_path = select_first([input_cram])

    call CramToBam {
      input:
        ref_fasta = select_first([ref_fasta]),
        ref_fasta_index = select_first([ref_fasta_index]),
        cram_file = select_first([input_cram]),
        output_basename = bam_from_cram_name,
        disk_size = ceil(cram_size * 6) + additional_disk
    }
  }

  File input_file = select_first([CramToBam.output_bam, input_bam])
  Float input_size = size(input_file, "GiB")

  if (!defined(output_map)) {
    call GenerateOutputMap {
      input:
        input_bam = input_file,
        unmapped_bam_suffix = unmapped_bam_suffix,
        disk_size = ceil(input_size) + additional_disk
    }
  }

  call SplitUpOutputMapFile {
    input:
      read_group_map_file = select_first([output_map, GenerateOutputMap.output_map])
  }

  scatter (rg_map_file in SplitUpOutputMapFile.rg_to_ubam_file) {
    call SplitOutUbamByReadGroup {
      input:
        input_bam = input_file,
        rg_to_ubam_file = rg_map_file,
        disk_size = ceil(input_size * 2) + additional_disk
    }

    String unmapped_bam_filename = basename(SplitOutUbamByReadGroup.output_bam)

    call RevertSam {
      input:
        input_bam = SplitOutUbamByReadGroup.output_bam,
        output_bam_filename = unmapped_bam_filename,
        disk_size = ceil(input_size * 3) + additional_disk
    }

    call SortSam {
      input:
        input_bam = RevertSam.output_bam,
        output_bam_filename = unmapped_bam_filename
    }

    Float unmapped_bam_size = size(SortSam.output_bam, "GiB")

    call ValidateSamFile {
      input:
        input_bam = SortSam.output_bam,
        report_filename = unmapped_bam_filename + ".validation_report",
        disk_size = ceil(unmapped_bam_size) + additional_disk
    }
  }

  output {
    Array[File] validation_report = ValidateSamFile.report
    Array[File] unmapped_bams = SortSam.output_bam
  }
  meta {
    allowNestedInputs: true
  }
}

task RevertSam {
  input {
    File input_bam
    String output_bam_filename
    Int disk_size
    Int memory_in_MiB = 3000
  }

  Int java_mem = memory_in_MiB - 1000
  Int max_heap = memory_in_MiB - 500

  command <<<
    java -Xms~{java_mem}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
    RevertSam \
    --INPUT ~{input_bam} \
    --OUTPUT ~{output_bam_filename} \
    --VALIDATION_STRINGENCY LENIENT \
    --ATTRIBUTE_TO_CLEAR FT \
    --ATTRIBUTE_TO_CLEAR CO \
    --ATTRIBUTE_TO_CLEAR PA \
    --ATTRIBUTE_TO_CLEAR OA \
    --ATTRIBUTE_TO_CLEAR XA \
    --SORT_ORDER coordinate

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory_in_MiB} MiB"
    preemptible: 3
  }

  output {
    File output_bam = output_bam_filename
  }
}

# This task is slower than converting straight from cram to bam (note we stream through sam format in between cram and bam)
# This is currently necessary due to a difference in the way the NM tag is calculated in order to produce a valid bam.
task CramToBam {
  input {
    File ref_fasta
    File ref_fasta_index
    File cram_file
    String output_basename
    Int disk_size
    Int memory_in_MiB = 7000
  }

  command <<<

    set -e
    set -o pipefail

    samtools view -h -T ~{ref_fasta} ~{cram_file} |
    samtools view -b -o ~{output_basename}.bam -
    samtools index -b ~{output_basename}.bam
    mv ~{output_basename}.bam.bai ~{output_basename}.bai

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
    cpu: 3
    memory: "~{memory_in_MiB} MiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }

  output {
    File output_bam = "~{output_basename}.bam"
    File output_bam_index = "~{output_basename}.bai"
  }
}

task GenerateOutputMap {
  input {
    File input_bam
    String unmapped_bam_suffix
    Int disk_size
    Int memory_in_MiB = 3000
  }

  command <<<

    set -e

    samtools view -H ~{input_bam} | grep '^@RG' | cut -f2 | sed s/ID:// > readgroups.txt

    echo -e "#READ_GROUP_ID\tOUTPUT" > output_map.tsv

    for rg in `cat readgroups.txt`; do
      echo -e "$rg\t$rg~{unmapped_bam_suffix}" >> output_map.tsv
    done

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory_in_MiB} MiB"
    preemptible: 3
  }

  output {
    File output_map = "output_map.tsv"
  }
}

task SplitUpOutputMapFile {
  input {
    File read_group_map_file
    Int disk_size = 10
    Int memory_in_MiB = 3000
  }

  command <<<
    mkdir rgtemp
    cd rgtemp

    # splits list of mappings into single files.  One line each.
    grep -v '^#' ~{read_group_map_file} | split -l 1 - rg_to_ubam_
  >>>

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory_in_MiB} MiB"
  }

  output {
    Array[File] rg_to_ubam_file = glob("rgtemp/rg_to_ubam_*")
  }
}

task SplitOutUbamByReadGroup {

  input {
    File input_bam
    File rg_to_ubam_file
    Int disk_size
    Int memory_in_MiB = 30000
  }

  Array[Array[String]] tmp = read_tsv(rg_to_ubam_file)

  command <<<
    echo "Read Group ~{tmp[0][0]} from ~{input_bam} is being written to ~{tmp[0][1]}"
    samtools view -b -h -r ~{tmp[0][0]} -o ~{tmp[0][1]} ~{input_bam}
  >>>

  output {
    File output_bam = tmp[0][1]
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
    cpu: 2
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory_in_MiB} MiB"
    preemptible: 3
  }
}

task ValidateSamFile {
  input {
    File input_bam
    String report_filename
    Int disk_size
    Int memory_in_MiB = 3000
  }

  Int java_mem = memory_in_MiB - 1000
  Int max_heap = memory_in_MiB - 500

  command <<<

    java -Xms~{java_mem}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
      ValidateSamFile \
      --INPUT ~{input_bam} \
      --OUTPUT ~{report_filename} \
      --MODE VERBOSE \
      --IS_BISULFITE_SEQUENCED false

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory_in_MiB} MiB"
    preemptible: 3
  }

  output {
    File report = "~{report_filename}"
  }
}

task SortSam {
  input {
    File input_bam
    String output_bam_filename
    Int memory_in_MiB = 7000
    Float sort_sam_disk_multiplier = 6
  }
  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20
  Int java_mem = memory_in_MiB - 1000
  Int max_heap = memory_in_MiB - 500

  command <<<
    java -Xms~{java_mem}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
    SortSam \
    --INPUT ~{input_bam} \
    --OUTPUT ~{output_bam_filename} \
    --SORT_ORDER queryname \
    --MAX_RECORDS_IN_RAM 300000

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory_in_MiB} MiB"
    preemptible: 3
  }

  output {
    File output_bam = output_bam_filename
  }
}
