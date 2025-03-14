version 1.0

# scATAC is now deprecated 2025-03-06

workflow scATAC {
    meta {
      description: "Processing of single-cell ATAC-seq data with the scATAC pipeline."
    }

    input {
        File input_fastq1
        File input_fastq2
        String input_id
        String genome_name
        File input_reference
        String output_bam = input_id + "_aligned.bam"
        String bin_size_list = "10000"
    }

    String pipeline_version = "1.3.3"

    parameter_meta {
        input_fastq1: "read 1 input fastq, the read names must be tagged with the cellular barcodes"
        input_fastq2: "read 2 input fastq, the read names must be tagged with the cellular barcodes"
        input_id: "name of the sample, used to name the outputs"
        input_reference: "tar file with BWA reference, generated with the build_bwa_reference pipeline"
        output_bam: "output BAM file name"
        genome_name: "name of the genome for scATAC"
        bin_size_list: "space separated list of bins to generate"
    }

    call AlignPairedEnd {
        input:
            input_fastq1 = input_fastq1,
            input_fastq2 = input_fastq2,
            input_reference = input_reference,
            output_bam = output_bam
    }

    call SnapPre {
        input:
            input_bam = AlignPairedEnd.aligned_bam,
            output_snap_basename = input_id + '.snap',
            genome_name = genome_name,
            input_reference = input_reference,
    }

    call SnapCellByBin {
        input:
            snap_input = SnapPre.output_snap,
            bin_size_list = bin_size_list,
            snap_output_name = input_id + '.snap'
    }

    call MakeCompliantBAM {
        input:
            input_bam = AlignPairedEnd.aligned_bam,
            output_bam_filename = input_id + '.bam'
    }

    call BreakoutSnap {
        input:
            snap_input = SnapCellByBin.output_snap,
            bin_size_list = bin_size_list,
            input_id = input_id
    }

    output {
        File output_snap_qc = SnapPre.output_snap_qc
        File output_snap = SnapCellByBin.output_snap
        File output_aligned_bam = MakeCompliantBAM.output_bam
        File breakout_barcodes = BreakoutSnap.barcodes
        File breakout_fragments = BreakoutSnap.fragments
        File breakout_binCoordinates = BreakoutSnap.binCoordinates
        File breakout_binCounts = BreakoutSnap.binCounts
        File breakout_barcodesSection = BreakoutSnap.barcodesSection
        String output_pipeline_version = pipeline_version
    }

}

task AlignPairedEnd {
    input {
        File input_fastq1
        File input_fastq2
        File input_reference
        String reference_unpack_name = "genome/genome.fa"
        String output_bam
        Int min_cov = 0
        String docker_image = "us.gcr.io/broad-gotc-prod/snaptools-bwa:1.0.0-1.4.8-0.7.17-1690310027"
        Int machine_mem_mb = 16000
        Int cpu = 16
        Int disk = ceil(2*(size(input_fastq1, "GiB") + size(input_fastq2, "GiB") + size(input_reference, "GiB"))) + 100
        Int preemptible = 3
    }

    parameter_meta {
      input_fastq1: "read 1 fastq file"
      input_fastq2: "read 2 fastq file"
      input_reference: "input reference bundle"
      reference_unpack_name: "name of the reference input file after decompression, used only if default is changed"
      output_bam: "name of output bam file"
      min_cov: "--min-cov parameter for snaptools align-paired-end (default: 0)"
    }


    command {
        set -euo pipefail

        # Make temp directory
        declare -r TEMP_DIR=`mktemp -d tmpdir_XXXXXX`

        # Unpack the reference
        tar xf ~{input_reference}

        # Run snaptools alignment
        snaptools align-paired-end \
            --input-reference=~{reference_unpack_name} \
            --input-fastq1=~{input_fastq1} \
            --input-fastq2=~{input_fastq2} \
            --output-bam=~{output_bam} \
            --aligner=bwa \
            --path-to-aligner=/usr/local/bin/ \
            --read-fastq-command=zcat \
            --min-cov=~{min_cov} \
            --num-threads=~{cpu} \
            --tmp-folder=$TEMP_DIR \
            --overwrite=TRUE \
            --if-sort=True
    }

    output {
        File aligned_bam = output_bam
    }

    runtime {
        docker: docker_image
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        disk: disk + " GB" # TES
        cpu: cpu
        preemptible: preemptible
    }
}

task SnapPre {
    input {
        File input_bam
        String output_snap_basename
        String genome_name
        String genome_size_file = "genome/chrom.sizes"
        String docker_image = "us.gcr.io/broad-gotc-prod/snaptools-bwa:1.0.0-1.4.8-0.7.17-1690310027"
        File input_reference
        Int cpu = 1
        Int machine_mem_mb = 16000
        Int disk = 500
        Int preemptible = 3
    }

    parameter_meta {
       input_bam: "input bam file"
       output_snap_basename: "prefix name of the output bam file"
       genome_name: "name of the genome, currently not parsed but saved in output"
       genome_size_file: "name of the chrom.sizes file after unpacking of the genome reference"
       docker_image: "docker image used"
       input_reference: "input reference tar file"
    }



    command {
        set -euo pipefail

        tar xf ~{input_reference}

        # Does the main counting
        snaptools snap-pre \
            --input-file=~{input_bam} \
            --output-snap=~{output_snap_basename} \
            --genome-name=~{genome_name} \
            --genome-size=~{genome_size_file} \
            --min-mapq=30  \
            --min-flen=0  \
            --max-flen=1000  \
            --keep-chrm=TRUE  \
            --keep-single=TRUE  \
            --keep-secondary=False  \
            --overwrite=True  \
            --max-num=1000000  \
            --min-cov=100  \
            --verbose=True
    }

    output {
        File output_snap = output_snap_basename
        File output_snap_qc = output_snap_basename + ".qc"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        disk: disk + " GB" # TES
        cpu: cpu
        preemptible: preemptible
    }
}

task SnapCellByBin {
    input {
        File snap_input
        String bin_size_list
        String snap_output_name
        String docker_image = "us.gcr.io/broad-gotc-prod/snaptools-bwa:1.0.0-1.4.8-0.7.17-1690310027"
        Int cpu = 1
        Int machine_mem_mb = 16000
        Int disk = 500
        Int preemptible = 3
    }

    parameter_meta {
       snap_input: "input snap file to generate bins for"
       bin_size_list: "list of bin sizes to generate"
       snap_output_name: "name of the output snap file"
       docker_image: "docker image to use"
    }

    command {
        set -euo pipefail

        mv ~{snap_input} ~{snap_output_name}

        # This is mutating the file in-place
        snaptools snap-add-bmat  \
            --snap-file ~{snap_output_name}  \
            --bin-size-list ~{bin_size_list}  \
            --verbose=True
    }

    output {
        File output_snap = snap_output_name
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        disk: disk + " GB" # TES
        preemptible: preemptible
    }
}

task MakeCompliantBAM {
    input {
        File input_bam
        String output_bam_filename
        String docker_image = "us.gcr.io/broad-gotc-prod/warp-tools:1.0.1-1690997141"
        Int cpu = 1
        Int disk =  ceil(3 * (size(input_bam, "GiB"))) + 100
        Int machine_mem_mb = 4000
        Int preemptible = 3

    }

    parameter_meta {
        input_bam: "input bam file"
        output_bam_filename: "name of output bam file"
        docker_image: "docker image to use"
    }

    command {
        set -euo pipefail

        /warptools/scripts/makeCompliantBAM.py --input-bam ~{input_bam} --output-bam ~{output_bam_filename}
    }

    output {
        File output_bam = output_bam_filename
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        disk: disk + " GB" # TES
        preemptible: preemptible
    }
}

task BreakoutSnap {
    input {
        File snap_input
        String docker_image = "us.gcr.io/broad-gotc-prod/warp-tools:1.0.1-1690997141"
        String bin_size_list
        String input_id
        Int preemptible = 3
        Int disk =  ceil(10 * (if size(snap_input, "GiB") < 1 then 1 else size(snap_input, "GiB") )) + 100
        Int machine_mem_mb = 16000
        Int cpu = 1
    }

    parameter_meta {
        snap_input: "input snap file to use"
        docker_image: "docker image to use"
        bin_size_list: "space separated list of bins to generate"
        input_id : "name of the sample, used to name the outputs"
    }

    command {
        set -euo pipefail
        mkdir output
        python3 /warptools/scripts/breakoutSnap.py --input ~{snap_input} \
            --output-prefix output/~{input_id}_
    }

    output {
        File barcodes = 'output/~{input_id}_barcodes.csv'
        File fragments = 'output/~{input_id}_fragments.csv'
        File binCoordinates = 'output/~{input_id}_binCoordinates_~{bin_size_list}.csv'
        File binCounts = 'output/~{input_id}_binCounts_~{bin_size_list}.csv'
        File barcodesSection = 'output/~{input_id}_barcodesSection.csv'
    }

    runtime {
        docker: docker_image
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        disk: disk + " GB" # TES
        cpu: cpu
        preemptible: preemptible
    }
}
