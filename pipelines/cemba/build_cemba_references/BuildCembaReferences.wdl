version 1.0

workflow BuildCembaReferences {
  input {
    File reference_fasta
    File? monitoring_script
  }

  # prepare reference fasta for building indexes with bowtie2
  call BuildBisulfiteReferences as Convert {
    input:
      fasta_input = reference_fasta,
      monitoring_script = monitoring_script
  }

  # build with the forward reference fasta
  call Bowtie2Build as IndexForward {
    input:
      fasta_input = Convert.fwd_converted_reference_fasta_output,
      index_prefix = "BS_CT",
      monitoring_script = monitoring_script
  }

  # build with the reverse reference fasta
  call Bowtie2Build as IndexReverse {
    input:
      fasta_input = Convert.rev_converted_reference_fasta_output,
      index_prefix = "BS_GA",
      monitoring_script = monitoring_script
  }

  call CreateReferenceDictionary {
    input:
      reference_fasta = reference_fasta,
      monitoring_script = monitoring_script
  }

  call CreateReferenceFastaIndex {
    input:
      reference_fasta = reference_fasta,
      monitoring_script = monitoring_script
  }

  # save the outputs for using with bismark/bowtie2
  output {
    File reference_fasta_dict =  CreateReferenceDictionary.ref_dict_output
    File reference_fasta_index = CreateReferenceFastaIndex.ref_index_output
    File fwd_converted_reference_fasta = Convert.fwd_converted_reference_fasta_output
    File rev_converted_reference_fasta = Convert.rev_converted_reference_fasta_output
    Array[File] fwd_bowtie2_index_files = IndexForward.bowtie2_index_files
    Array[File] rev_bowtie2_index_files = IndexReverse.bowtie2_index_files
  }
}


# converts bases of reference fasta to build indexes of bowtie2
task BuildBisulfiteReferences {
  input {
    File fasta_input
    File? monitoring_script
  }

  # input file size
  Float input_size = size(fasta_input, "GB")

  # output names of converted fasta's
  String fwd_converted_reference_fasta_output_name = "genome_mfa.CT_conversion.fa"
  String rev_converted_reference_fasta_output_name = "genome_mfa.GA_conversion.fa"

  # custom script for reading in fasta and changing c's to t's (fwd) and g's to a's (rev)
  command <<<
    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    python /build_bisulfite_references.py \
      --input-fasta ~{fasta_input} \
      --forward-convert-out ~{fwd_converted_reference_fasta_output_name} \
      --reverse-convert-out ~{rev_converted_reference_fasta_output_name}
  >>>

  runtime {
    docker: "quay.io/broadinstitute/bisulfite-references:1.0"
    disks: "local-disk " + ceil(3.5 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3.5 GB"
  }

  output {
    File fwd_converted_reference_fasta_output = fwd_converted_reference_fasta_output_name
    File rev_converted_reference_fasta_output = rev_converted_reference_fasta_output_name
    File monitoring_log = "monitoring.log"
  }
}

# build with the converted fasta's
task Bowtie2Build {
  input {
    File fasta_input
    String index_prefix
    File? monitoring_script
  }

  # input file size
  Float input_size = size(fasta_input, "GB")

  # use bowtie2 to build indexes
  command <<<
    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    bowtie2-build \
      -f ~{fasta_input} ~{index_prefix}
  >>>

  runtime {
    docker: "quay.io/broadinstitute/bowtie2:2.3.4.3"
    disks: "local-disk " + ceil(3 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "7 GB"
  }
  
  output {
    Array[File] bowtie2_index_files = glob(index_prefix + "*")
    File monitoring_log = "monitoring.log"
  }
}

# create reference dictionary
task CreateReferenceDictionary {
  input {
    File reference_fasta
    File? monitoring_script
  }

  # input file size
  Float input_size = size(reference_fasta, "GB")

  # output names for read1 bam with extracted barcodes
  String ref_dict_output_name = basename(basename(reference_fasta, ".fa"), ".fasta") + ".dict"

  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    # create a reference dict
    java -jar /picard-tools/picard.jar CreateSequenceDictionary \
      REFERENCE=~{reference_fasta} \
      OUTPUT=~{ref_dict_output_name}

  >>>

  # use docker image for given tool cutadapat
  runtime {
    docker: "quay.io/broadinstitute/picard:2.18.23"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 2 * input file size
    disks: "local-disk " + ceil(2 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "4 GB"
  }

  output {
    File ref_dict_output = ref_dict_output_name
    File monitoring_log = "monitoring.log"
  }
}

# create reference dictionary
task CreateReferenceFastaIndex {
  input {
    File reference_fasta
    File? monitoring_script
  }

  # input file size
  Float input_size = size(reference_fasta, "GB")

  # output names for read1 bam with extracted barcodes
  String ref_index_output_name = basename(reference_fasta) + ".fai"

  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    declare -a FASTA_TEMP=$(mktemp fasta_XXXXXX)
    cp ~{reference_fasta} $FASTA_TEMP

    # create a reference index
    samtools faidx \
        $FASTA_TEMP

    mv $FASTA_TEMP.fai ~{ref_index_output_name}
  >>>

  # use docker image for given tool cutadapat
  runtime {
    docker: "quay.io/broadinstitute/samtools:1.9"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 2.25 * input file size
    disks: "local-disk " + ceil(2.25 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3.5 GB"
  }

  output {
    File ref_index_output = ref_index_output_name
    File monitoring_log = "monitoring.log"
  }
}
