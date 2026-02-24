version 1.0

task count {
  input {
    Array[String] fastq_paths
    Array[String] pucks
    String input_id
    Int mem_GiB = 64
    Int disk_GiB = 128
    Int nthreads = 1
    String docker
    String? billing_project
  }

  String gsutil_requester_pays_flag = if defined(billing_project) then "-u ${billing_project}" else ""
  String gcloud_requester_pays_flag = if defined(billing_project) then "--billing-project ${billing_project}" else ""

  command <<<
    set -euo pipefail
    set -x

    echo "<< starting spatial-count >>"

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Download the script -- put this script into a docker 
    # This needs to be changed when Matthew finalizes changes to these scripts
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/d89176cf21e072fe8b5aad3a1454ad194fca7c9a/slide-tags/spatial-count.jl

    echo "FASTQs: ~{length(fastq_paths)} paths provided"
    echo "Pucks: ~{length(pucks)} puck(s) provided"

    # Assert that the fastqs exist
    fastqs=(~{sep=' ' fastq_paths})
    for fastq in "${fastqs[@]}" ; do
        if ! gsutil ~{gsutil_requester_pays_flag} stat "$fastq" &> /dev/null ; then
            echo "ERROR: gsutil stat command failed on fastq $fastq"
            exit 1
        fi
    done

    # Download the fastqs
    echo "Downloading fastqs:"
    mkdir fastqs
    gcloud storage  ~{gcloud_requester_pays_flag} cp ~{sep=' ' fastq_paths} fastqs

    # Assert that the pucks exist
    pucks=(~{sep=' ' pucks})
    for puck in "${pucks[@]}" ; do
        if ! gsutil ~{gsutil_requester_pays_flag} stat "$puck" &> /dev/null ; then
            echo "ERROR: gsutil stat command failed on puck $puck"
            exit 1
        fi
    done

    # Download the pucks
    echo "Downloading pucks:"
    mkdir pucks
    gcloud storage  ~{gcloud_requester_pays_flag} cp ~{sep=' ' pucks} pucks

    # Run the script
    echo ; echo "Running spatial-count.jl"
    ## julia --threads=4 /spatial-count.jl fastqs pucks .
    julia --threads=4 spatial-count.jl fastqs pucks .
    
    if [[ -f SBcounts.h5 ]] ; then
        echo ; echo "Success, uploading counts"
        echo "true" > DONE
    else
        echo ; echo "ERROR: CANNOT FIND: SBcounts.h5"
    fi

    echo; echo "Writing logs:"
    echo; echo "fastqs size:"; du -sh fastqs
    echo; echo "pucks size:"; du -sh pucks
    echo; echo "output size:"; du -sh SBcounts.h5
    echo; echo "FREE SPACE:"; df -h

    mv SBcounts.h5 ~{input_id}_SBcounts.h5 
    cat stdout stderr > ~{input_id}_spatial-count.log
    echo "<< completed spatial-count >>"

  >>>
  
  output {
    File sb_counts = "~{input_id}_SBcounts.h5"
    File spatial_log = "~{input_id}_spatial-count.log"
  }
  
  runtime {
    docker: docker
    memory: "~{mem_GiB} GB"
    disks: "local-disk ~{disk_GiB} SSD"
    cpu: nthreads
  }
}

