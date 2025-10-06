version 1.0

task FastqProcessing {
  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    Array[File]? i1_fastq
    File whitelist
    Int chemistry
    String sample_id
    String read_struct

    #using the latest build of warp-tools in GCR
    String warp_tools_docker_path

    #runtime values
    Int machine_mem_mb = 40000
    Int cpu = 16   
    #TODO decided cpu
    # estimate that bam is approximately equal in size to fastq, add 20% buffer
    Int disk = ceil(size(r1_fastq, "GiB")*3 + size(r2_fastq, "GiB")*3) + 500

    Int preemptible = 3
  }

  meta {
    description: "Converts a set of fastq files to unaligned bam file, also corrects barcodes and partitions the alignments by barcodes."
  }

  parameter_meta {
    r1_fastq: "input fastq file"
    r2_fastq: "input fastq file"
    i1_fastq: "(optional) input fastq file"
    read_struct: "read structure for the 10x chemistry. This automatically selected in the checkInputs task"
    whitelist: "10x genomics cell barcode whitelist"
    chemistry: "chemistry employed, currently can be tenX_v2 or tenX_v3, the latter implies NO feature barcodes"
    sample_id: "name of sample matching this file, inserted into read group header"
    warp_tools_docker_path: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    FASTQS=$(python3 <<CODE
    def rename_file(filename):
        import shutil
        import gzip
        import re
         
        iscompressed = True
        with gzip.open(filename, 'rt') as fin:
           try:
               _ = fin.readline()
           except:
               iscompressed = False

        basename = re.sub(r'.gz$', '', filename)
        basename = re.sub(r'.fastq$', '', basename)
 
        if iscompressed:
            # if it is already compressed then add an extension .fastq.gz
            newname = basename + ".fastq.gz" 
        else: 
            # otherwise, add just the .fastq extension
            newname = basename + ".fastq"

        if filename != newname:
            # safe to rename since the old and the new names are different
            shutil.move(filename, newname)

        return newname
    optstring = ""
     
    r1_fastqs = [ "${sep='", "' r1_fastq}" ]
    r2_fastqs = [ "${sep='", "' r2_fastq}" ]
    i1_fastqs = [ "${sep='", "' i1_fastq}" ]
    for fastq in r1_fastqs:
        if fastq.strip(): 
            optstring += " --R1 " + rename_file(fastq)
    for fastq in r2_fastqs:
        if fastq.strip(): 
            optstring += " --R2 " + rename_file(fastq)
    for fastq in i1_fastqs:
        if fastq.strip(): 
            optstring += " --I1 " + rename_file(fastq)
    print(optstring)
    CODE)

    # use the right UMI length depending on the chemistry
    if [ "~{chemistry}" == "2" ]; then
        ## V2
        UMILENGTH=10
    elif [ "~{chemistry}" == "3" ]; then
        ## V3
        UMILENGTH=12
    else
        echo Error: unknown chemistry value: "~{chemistry}"
        exit 1;
    fi

    fastqprocess \
        --bam-size 30.0 \
        --sample-id "~{sample_id}" \
        $FASTQS \
        --white-list "~{whitelist}" \
        --read-structure "~{read_struct}" \
        --output-format FASTQ
  }
  
  runtime {
    docker: warp_tools_docker_path
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    Array[File] fastq_R1_output_array = glob("fastq_R1_*")
    Array[File] fastq_R2_output_array = glob("fastq_R2_*")
  }
}

task FastqProcessingSlidSeq {

  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    Array[File]? i1_fastq
    String read_structure
    String sample_id
    File whitelist


    # Runtime attributes
    String docker =  "us.gcr.io/broad-gotc-prod/warp-tools:2.6.1"
    Int cpu = 16
    Int machine_mb = 40000
    Int disk = ceil(size(r1_fastq, "GiB")*3 + size(r2_fastq, "GiB")*3) + 50
    Int preemptible = 3
  }

  meta {
    description: "Converts a set of fastq files to unaligned bam file, also corrects barcodes and partitions the alignments by barcodes. Allows for variable barcode and umi lengths as input"
  }

  parameter_meta {
        r1_fastq: "Array of Read 1 FASTQ files - forward read, contains cell barcodes and molecule barcodes"
        r2_fastq: "Array of Read 2 FASTQ files - reverse read, contains cDNA fragment generated from captured mRNA"
        i1_fastq: "(optional) Array of i1 FASTQ files - index read, for demultiplexing of multiple samples on one flow cell."
        sample_id: "Name of sample matching this file, inserted into read group header"
        read_structure: "A string that specifies UMI (M) and Barcode (C) positions in the Read 1 fastq"

  }

  command {

    command {
    set -e

    FASTQS=$(python3 <<CODE
    def rename_file(filename):
        import shutil
        import gzip
        import re

        iscompressed = True
        with gzip.open(filename, 'rt') as fin:
          try:
              _ = fin.readline()
          except:
              iscompressed = False

        basename = re.sub(r'.gz$', '', filename)
        basename = re.sub(r'.fastq$', '', basename)

        if iscompressed:
            # if it is already compressed then add an extension .fastq.gz
            newname = basename + ".fastq.gz"
        else:
            # otherwis, add just the .fastq extension
            newname = basename + ".fastq"

        if filename != newname:
            # safe to rename since the old and the new names are different
            shutil.move(filename, newname)

        return newname
    optstring = ""

    r1_fastqs = [ "${sep='", "' r1_fastq}" ]
    r2_fastqs = [ "${sep='", "' r2_fastq}" ]
    i1_fastqs = [ "${sep='", "' i1_fastq}" ]
    for fastq in r1_fastqs:
        if fastq.strip():
            optstring += " --R1 " + rename_file(fastq)
    for fastq in r2_fastqs:
        if fastq.strip():
            optstring += " --R2 " + rename_file(fastq)
    for fastq in i1_fastqs:
        if fastq.strip():
            optstring += " --I1 " + rename_file(fastq)
    print(optstring)
    CODE)
    cut -f 1 ~{whitelist} > WhiteList.txt


    fastq_slideseq  \
      --bam-size 30.0 \
      --white-list WhiteList.txt \
      --read-structure "~{read_structure}" \
      --sample-id "~{sample_id}" \
      --output-format FASTQ \
      $FASTQS
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "${machine_mb} MiB"
    disks: "local-disk ${disk} HDD"
    preemptible: preemptible
  }
  
  output {
    Array[File] fastq_R1_output_array = glob("fastq_R1_*")
    Array[File] fastq_R2_output_array = glob("fastq_R2_*")
  }
}

task FastqProcessATAC {

    input {
        Array[File] read1_fastq
        Array[File] read3_fastq
        Array[File] barcodes_fastq
        String read_structure = "16C"
        String barcode_orientation = "FIRST_BP_RC"
        String output_base_name
        File whitelist
        String barcode_index1 = basename(barcodes_fastq[0])
        String docker_path

        # Runtime attributes [?]
        Int mem_size = 5
        Int cpu = 16
        # TODO decided cpu
        # estimate that bam is approximately equal in size to fastq, add 20% buffer
        Int disk_size = ceil(2 * ( size(read1_fastq, "GiB") + size(read3_fastq, "GiB") + size(barcodes_fastq, "GiB") )) + 400
        Int preemptible = 3

        # Additional parameters for fastqprocess
        Int num_output_files
    }

    meta {
        description: "Converts a set of fastq files to unaligned bam file/fastq file, also corrects barcodes and partitions the alignments by barcodes. Allows for variable barcode and umi lengths as input, if applicable."
    }

    parameter_meta {
        read1_fastq: "Array of read 1 FASTQ files of paired reads -- forward reads"
        read3_fastq: "Array of read 3 FASTQ files of paired reads -- reverse reads"
        barcodes_fastq: "Array of read 2 FASTQ files which contains the cellular barcodes"
        output_base_name: "Name of sample matching this file, inserted into read group header"
        read_structure: "A string that specifies the barcode (C) positions in the Read 2 fastq"
        barcode_orientation: "A string that specifies the orientation of barcode needed for scATAC data. The default is FIRST_BP. Other options include LAST_BP, FIRST_BP_RC or LAST_BP_RC."
        whitelist: "10x genomics cell barcode whitelist for scATAC"
        docker_path: "The docker image path containing the runtime environment for this task"
        mem_size: "(optional) the amount of memory (MiB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
        disk_size: "(optional) the amount of disk space (GiB) to provision for this task"
        num_output_files: "(optional) the number of output fastq file shards to produce. if this is set to > 0, bam_size is ignored."
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    }

    command <<<

        set -e
        echo "Num of output files"
        echo ~{num_output_files} 
        
        declare -a FASTQ1_ARRAY=(~{sep=' ' read1_fastq})
        declare -a FASTQ2_ARRAY=(~{sep=' ' barcodes_fastq})
        declare -a FASTQ3_ARRAY=(~{sep=' ' read3_fastq})

        read1_fastq_files=`printf '%s ' "${FASTQ1_ARRAY[@]}"; echo`
        read2_fastq_files=`printf '%s ' "${FASTQ2_ARRAY[@]}"; echo`
        read3_fastq_files=`printf '%s ' "${FASTQ3_ARRAY[@]}"; echo`

        echo $read1_fastq_files
        # Make downsample fq for barcode orientation check of R2 barcodes
        mkdir -p input_fastqs

        # Function to move files into the input_fastqs directory
        move_files_to_input_dir() {
            local -n array=$1  # Reference to the array passed as argument
            local destination_dir=$2

            for file in "${array[@]}"; do
                if [ -f "$file" ]; then  # Check if file exists
                    echo "Moving $file to $destination_dir"
                    mv "$file" "$destination_dir"
                else
                    echo "File $file not found"
                fi
            done
        }

        # Move files from FASTQ1_ARRAY to input_fastqs directory
        move_files_to_input_dir FASTQ1_ARRAY input_fastqs

        # Move files from FASTQ2_ARRAY to input_fastqs directory
        move_files_to_input_dir FASTQ2_ARRAY input_fastqs

        # Move files from FASTQ3_ARRAY to input_fastqs directory
        move_files_to_input_dir FASTQ3_ARRAY input_fastqs

        echo "All files moved to input_fastqs directory"

        #gcloud storage cp $read1_fastq_files /cromwell_root/input_fastqs
        #gcloud storage cp $read2_fastq_files /cromwell_root/input_fastqs
        #gcloud storage cp $read3_fastq_files /cromwell_root/input_fastqs

        path="input_fastqs/"
        barcode_index="~{barcode_index1}"
        file="${path}${barcode_index}"
        zcat "$file" | sed -n '2~4p' | shuf -n 1000 > downsample.fq
        head -n 1 downsample.fq
        # barcodes R2
        R1_FILES_CONCAT=""
        for fastq in "${FASTQ2_ARRAY[@]}"
        do
            BASE=`basename $fastq`
            BASE=`echo --R1 input_fastqs/$BASE`
            R1_FILES_CONCAT+="$BASE "
        done
        echo $R1_FILES_CONCAT

        # R1
        R2_FILES_CONCAT=""
        for fastq in "${FASTQ1_ARRAY[@]}"
        do
            BASE=`basename $fastq`
            BASE=`echo --R2 input_fastqs/$BASE`
            R2_FILES_CONCAT+="$BASE "
        done
        echo $R2_FILES_CONCAT
        
        # R3
        R3_FILES_CONCAT=""
        for fastq in "${FASTQ3_ARRAY[@]}"
        do
            BASE=`basename $fastq`
            BASE=`echo --R3 input_fastqs/$BASE`
            R3_FILES_CONCAT+="$BASE "
        done
        echo $R3_FILES_CONCAT

        python3 /warptools/scripts/dynamic-barcode-orientation.py downsample.fq "~{whitelist}" best_match.txt
        
        cat best_match.txt
        barcode_choice=$(<best_match.txt)
        echo $barcode_choice

        # Call fastq process
        # outputs fastq files where the corrected barcode is in the read name

        fastqprocess \
        --num-output-files ~{num_output_files} \
        --sample-id "~{output_base_name}" \
        $R1_FILES_CONCAT \
        $R2_FILES_CONCAT \
        $R3_FILES_CONCAT \
        --white-list "~{whitelist}" \
        --output-format "FASTQ" \
        --barcode-orientation $barcode_choice \
        --read-structure "~{read_structure}"

    >>>

    runtime {
        docker: docker_path
        cpu: cpu
        memory: "${mem_size} MiB"
        disks: "local-disk ${disk_size} HDD"
        preemptible: preemptible
    }

    output {
        Array[File] fastq_R1_output_array = glob("fastq_R1_*")
        Array[File] fastq_R3_output_array = glob("fastq_R3_*")
    }
}

