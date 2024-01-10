version 1.0
task PairedTagDemultiplex {
    input {
        File read1_fastq
        File read3_fastq
        File barcodes_fastq
        String input_id
        Boolean preindex
        File whitelist
        String docker = "us.gcr.io/broad-gotc-prod/upstools:1.2.0-2023.03.03-1704723060"
        Int cpu = 1
        Int disk_size = ceil(2 * (size(read1_fastq, "GiB") + size(read3_fastq, "GiB") + size(barcodes_fastq, "GiB") )) + 400
        Int preemptible = 3
        Int mem_size = 8
    }
    meta {
        description: "Checks read2 FASTQ length and orientation and performs trimming."
    }
    parameter_meta {
        read1_fastq: "read 1 FASTQ files of paired reads -- forward reads"
        read3_fastq: "read 3 FASTQ files of paired reads -- reverse reads"
        barcodes_fastq: "read 2 FASTQ files which contains the cellular barcodes"
        preindex: "Boolean for whether data has a sample barcode that needs to be demultiplexed"
        whitelist: "Atac whitelist for 10x multiome data"
        input_id: "Input ID to demarcate sample"
        docker: "(optional) the docker image containing the runtime environment for this task"
        mem_size: "(optional) the amount of memory (MiB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
        disk_size: "(optional) the amount of disk space (GiB) to provision for this task"
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"        
    }
    command <<<
        set -e
        ## Need to gunzip the r1_fastq
        zcat ~{barcodes_fastq} | head -n2 > r2.fastq
        FASTQ=r2.fastq
        echo 'this is the fastq:' $FASTQ
        R2=$(awk 'NR==2' $FASTQ)
        COUNT=$(echo ${#R2})
        echo 'this is the read:' $R2
        echo 'this is the barcode count:' $COUNT
        echo "Renaming files for UPS tools"
        mv ~{read1_fastq} "~{input_id}_R1.fq.gz"
        mv ~{barcodes_fastq} "~{input_id}_R2.fq.gz"
        mv ~{read3_fastq} "~{input_id}_R3.fq.gz"
        echo performing read2 length and orientation checks 
        if [[ $COUNT == 27 && ~{preindex} == "false" ]]
          then
          pass="true"
          echo "Preindex is false and length is 27 bp"
          echo "Trimming first 3 bp with UPStools"
          upstools trimfq ~{input_id}_R2.fq.gz 4 26
          echo "Running orientation check"
          file="~{input_id}_R2_trim.fq.gz"
          zcat "$file" | sed -n '2~4p' | shuf -n 1000 > downsample.fq
          head -n 1 downsample.fq
          python3 /upstools/pyscripts/dynamic-barcode-orientation.py downsample.fq ~{whitelist} best_match.txt
          cat best_match.txt
          barcode_choice=$(<best_match.txt)
          echo "Barcode choice is: "
          echo $barcode_choice
          if [[ $barcode_choice == "FIRST_BP_RC" ]]; then
            pass="true"
          else
            pass="false"
            echo "Incorrect barcode orientation"
          fi
          mv "~{input_id}_R2_trim.fq.gz" "~{input_id}_R2.fq.gz"

        elif [[ $COUNT == 27 && ~{preindex} == "true" ]]
          then
          echo "Count is 27 bp because of preindex"
          echo "Running demultiplexing with UPStools"
          upstools sepType_DPT ~{input_id} 3
          echo "Running orientation check"
          file="~{input_id}_R2_prefix.fq.gz"
          zcat "$file" | sed -n '2~4p' | shuf -n 1000 > downsample.fq
          head -n 1 downsample.fq
          python3 /upstools/pyscripts/dynamic-barcode-orientation.py downsample.fq ~{whitelist} best_match.txt
          cat best_match.txt
          barcode_choice=$(<best_match.txt)
          echo "Barcode choice is: "
          echo $barcode_choice
          if [[ $barcode_choice == "FIRST_BP_RC" ]]; then
            pass="true"
          else
            pass="false"
            echo "Incorrect barcode orientation"
          fi
          # rename files to original name
          mv "~{input_id}_R2_prefix.fq.gz" "~{input_id}_R2.fq.gz"
          mv "~{input_id}_R1_prefix.fq.gz" "~{input_id}_R1.fq.gz"
          mv "~{input_id}_R3_prefix.fq.gz" "~{input_id}_R3.fq.gz"
        elif [[ $COUNT == 24 && ~{preindex} == "false" ]]
          then
          pass="true"
        echo "FASTQ has correct index length, no modification necessary"
        else
          echo "Length of read2 is not expected length; ending pipeline run"
          pass="false"
        fi
        if [[ $pass == "true" ]]
          then
          exit 0;
        else
          exit 1;
        fi
        exit 0;      
    >>>
    
    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_size} GiB"
        disks: "local-disk ${disk_size} HDD"
        preemptible: preemptible        
    }

    output {
        File fastq1 = "~{input_id}_R1.fq.gz"
        File barcodes = "~{input_id}_R2.fq.gz"
        File fastq3 = "~{input_id}_R3.fq.gz"
    }
}

task AddBBTag {
    input {
        File bam
        String input_id

        # using the latest build of upstools docker in GCR
        String docker = "us.gcr.io/broad-gotc-prod/upstools:1.0.0-2023.03.03-1704300311"

        # Runtime attributes
        Int mem_size = 8
        Int cpu = 1
        # TODO decided cpu
        # estimate that bam is approximately equal in size to fastq, add 20% buffer
        Int disk_size = ceil(2 * ( size(bam, "GiB"))) + 100
        Int preemptible = 3
    }

    meta {
        description: "Demultiplexes paired-tag ATAC fastq files that have a 3 bp preindex and adds the index to readnames."
    }

    parameter_meta {
        bam: "BAM with aligned reads and barcode in the CB tag"
        input_id: "input ID"
        docker: "(optional) the docker image containing the runtime environment for this task"
        mem_size: "(optional) the amount of memory (MiB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
        disk_size: "(optional) the amount of disk space (GiB) to provision for this task"
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    }

    command <<<

        set -e
        echo "BAM file name is:"
        echo ~{bam}
        echo moving BAM
        mv ~{bam} ./~{input_id}.bam
        echo Running UPStools
        python3 /upstools/pyscripts/scifi.preindex_CB_to_BB.py --in ~{input_id}.bam
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_size} GiB"
        disks: "local-disk ${disk_size} HDD"
        preemptible: preemptible
    }

    output {
        File bb_bam = "~{input_id}.bam.BB.bam"
    }
}