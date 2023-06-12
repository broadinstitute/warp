version 1.0


workflow ATAC {
    meta {
        description: "Processing for single-cell ATAC-seq data from the level of raw fastq reads. This is the first step of the multiome pipeline. ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) is a technique used in molecular biology to assess genome-wide chromatin accessibility. This pipeline processes 10x Genomics Multiome ATAC FASTQ files."
        allowNestedInputs: true
    }

    input {
    File read1_fastq = "gs://broad-gotc-test-storage/nikelle_testing/small_test.fastq"
    File read3_fastq = "gs://broad-gotc-test-storage/nikelle_testing/small_test2.fastq"
    File tar_bwa_reference = "gs://fc-dd55e131-ef49-4d02-aa2a-20640daaae1e/submissions/8f0dd71a-b42f-4503-b839-3f146941758a/IndexRef/53a91851-1f6c-4ab9-af66-b338ffb28b5a/call-BwaMem2Index/GRCh38.primary_assembly.genome.bwamem2.fa.tar"
    String output_base_name = "scATAC"
    File monitoring_script = "gs://fc-51792410-8543-49ba-ad3f-9e274900879f/cromwell_monitoring_script2.sh"
    }

    String pipeline_version = "1.0.1"


        call BWAPairedEndAlignment {
            input:
                read1_fastq = read1_fastq,
                read3_fastq = read3_fastq,
                tar_bwa_reference = tar_bwa_reference,
                output_base_name = output_base_name,
                monitoring_script = monitoring_script
        }


    output {
        File bam_aligned_output = BWAPairedEndAlignment.bam_aligned_output

    }
}
# align the two trimmed fastq as paired end data using BWA
task BWAPairedEndAlignment {
  input {
    File read1_fastq
    File read3_fastq
    File tar_bwa_reference
    String read_group_id = "RG1"
    String read_group_sample_name = "RGSN1"
    String output_base_name
    String docker_image = "us.gcr.io/broad-gotc-prod/samtools-bwa-mem-2:1.0.0-2.2.1_x64-linux-1685469504"

    # script for monitoring tasks
    File monitoring_script

    # Runtime attributes
    Int disk_size = ceil(3.25 * (size(read1_fastq, "GiB") + size(read3_fastq, "GiB") + size(tar_bwa_reference, "GiB"))) + 200
    Int nthreads = 16
    Int mem_size = 40
  }

  parameter_meta {
    read1_fastq: "the trimmed read 1 fastq file containing sequencing reads as input for the aligner"
    read3_fastq: "the trimmed read 1 fastq file containing sequencing reads as input for the aligner"
    tar_bwa_reference: "the pre built tar file containing the reference fasta and cooresponding reference files for the BWA aligner"
    read_group_id: "the read group id to be added upon alignment"
    read_group_sample_name: "the read group sample to be added upon alignment"
    nthreads: "the number of threads to use during bwa alignment"
    mem_size: "the size of memory used during alignment"
    disk_size : "disk size used in bwa alignment step"
    output_base_name: "basename to be used for the output of the task"
    docker_image: "the docker image using BWA to be used (default: us.gcr.io/broad-gotc-prod/samtools-bwa-mem-2:1.0.0-2.2.1_x64-linux-1685469504)"
    monitoring_script : "script to monitor resource comsumption of tasks"
  }

  String bam_aligned_output_name = output_base_name + ".aligned.bam"

  # bwa and call samtools to convert sam to bam
  command <<<

    set -euo pipefail

    if [ ! -z "~{monitoring_script}" ]; then
    chmod a+x ~{monitoring_script}
    ~{monitoring_script} > monitoring.log &
    else
    echo "No monitoring script given as input" > monitoring.log &
    fi

    # prepare reference
    declare -r REF_DIR=$(mktemp -d genome_referenceXXXXXX)
    tar -xf "~{tar_bwa_reference}" -C $REF_DIR --strip-components 1
    rm "~{tar_bwa_reference}"

    # align w/ BWA: -t for number of cores
    bwa-mem2 \
    mem \
    -R "@RG\tID:RG1\tSM:RGSN1" \
    -t 16 \
    $REF_DIR/genome.fa \
    ~{read1_fastq} ~{read3_fastq} \
    | samtools view -h \
    | awk '{if ($0 ~ /^@/) {print $0} else {cr=substr($1, index($1, "CR:")+3); n=split($1, a,":CB:"); if (n == 2) {cb="CB:Z:"a[1]"\t";} else {cb="";} print($0 "\tCR:Z:" cr "\t" cb "XC:Z:" substr($1, index($1, "CB:")+3));}}' \
    | samtools view -b -o ~{bam_aligned_output_name}
  >>>

  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: nthreads
    memory: "${mem_size} GiB"
  }

  output {
    File bam_aligned_output = bam_aligned_output_name
    File monitoring_log = "monitoring.log"
  }
}

