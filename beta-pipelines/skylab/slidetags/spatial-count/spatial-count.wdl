version 1.0

task count {
  input {
    Array[String] fastq_paths
    Array[String] pucks
    Int mem_GiB
    Int disk_GiB
    String docker
  }
  command <<<
    set -euo pipefail
    set -x
    # Taken from https://github.com/MacoskoLab/Macosko-Pipelines/blob/main/spatial-count/spatial-count.wdl
    # Modified to include the outputs
    # Last commit: eaebb9060fb05ececd980fd62438487d07990596 
  
    echo "<< starting spatial-count >>"
    dstat --time --cpu --mem --disk --io --freespace --output spatial-count.usage &> /dev/null &

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Download the script -- put this script into a docker
    ### wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/main/spatial-count/spatial-count.jl

    echo "FASTQs: ~{length(fastq_paths)} paths provided"
    echo "Pucks: ~{length(pucks)} puck(s) provided"
    echo "Output directory: $count_output_path" ; echo

    # Assert that the fastqs exist
    fastqs=(~{sep=' ' fastq_paths})
    for fastq in "${fastqs[@]}" ; do
        if ! gsutil stat "$fastq" &> /dev/null ; then
            echo "ERROR: gsutil stat command failed on fastq $fastq"
            exit 1
        fi
    done

    # Download the fastqs
    echo "Downloading fastqs:"
    mkdir fastqs
    gcloud storage cp ~{sep=' ' fastq_paths} fastqs

    # Assert that the pucks exist
    pucks=(~{sep=' ' pucks})
    for puck in "${pucks[@]}" ; do
        if ! gsutil stat "$puck" &> /dev/null ; then
            echo "ERROR: gsutil stat command failed on puck $puck"
            exit 1
        fi
    done

    # Download the pucks
    echo "Downloading pucks:"
    mkdir pucks
    gcloud storage cp ~{sep=' ' pucks} pucks

    # Run the script
    echo ; echo "Running spatial-count.jl"
    julia --threads=4 /spatial-count.jl fastqs pucks

    if [[ -f SBcounts.h5 ]] ; then
        echo ; echo "Success, uploading counts"
        ####gcloud storage cp -r SBcounts.h5 "$count_output_path/SBcounts.h5"
        echo "true" > DONE
    else
        echo ; echo "ERROR: CANNOT FIND: SBcounts.h5"
    fi

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "fastqs size:"; du -sh fastqs
    echo; echo "pucks size:"; du -sh pucks
    echo; echo "output size:"; du -sh SBcounts.h5
    echo; echo "FREE SPACE:"; df -h
     
    echo "zipping logs"
    tar -xvf /cromwell_root/stdout /cromwell_root/stderr /cromwell_root/spatial-count.usage.tar.gz

    echo "<< completed spatial-count >>"
  >>>
  
  output {
    Boolean DONE = read_boolean("DONE")
    File sb_counts = "SBcounts.h5"
    File spatial_logs = "spatial-count.usage.tar.gz"

  }
  runtime {
    docker: docker
    memory: "~{mem_GiB} GB"
    disks: "local-disk ~{disk_GiB} SSD"
    cpu: 1
    preemptible: 0
  }
}

