version 1.0

task count {
  input {
    Array[String] rna_paths
    String sb_path
    Int mem_GiB  = 128
    Int disk_GiB = 128
    String docker
  }
  command <<<
    echo "<< starting spatial-count >>"
    dstat --time --cpu --mem --disk --io --freespace --output positioning.usage &> /dev/null &

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Download the scripts
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/positioning/run-positioning.R
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/positioning/load_matrix.R
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/positioning/positioning.R

    echo "RNA: ~{sep=' ' rna_paths}"
    echo "SB: ~{sb_path}"

    # Assert that the RNA files exist
    rnas=(~{sep=' ' rna_paths})
    for rna in "${rnas[@]}" ; do
        if ! gsutil stat "$rna" &> /dev/null ; then
            echo "ERROR: gsutil stat command failed on file $rna"
            exit 1
        fi
    done

    # Download the RNA
    echo "Downloading RNA:"
    mkdir RNA
    gcloud storage cp ~{sep=' ' rna_paths} RNA

    # Assert that the SB file exists
    if ! gsutil stat "~{sb_path}" &> /dev/null ; then
        echo "ERROR: gsutil stat command failed on file ~{sb_path}"
        exit 1
    fi

    # Download the SB
    echo "Downloading SB:"
    mkdir SB
    gcloud storage cp ~{sb_path} SB

    # Run the script
    echo ; echo "Running run-positioning.R"
    Rscript run-positioning.R RNA SB output

    # Upload the results
    gcloud storage cp -r output/* "$count_output_path"

    if [[ -f output/seurat.qs ]] ; then
        echo "true" > DONE
    else
        echo ; echo "ERROR: CANNOT FIND: seurat.qs"
    fi

    echo; echo "Writing logs:"
    kill $(ps aux | fgrep dstat | fgrep -v grep | awk '{print $2}')
    echo; echo "RNA size:"; du -sh RNA
    echo; echo "SB size:"; du -sh SB
    echo; echo "output size:"; du -sh output
    echo; echo "FREE SPACE:"; df -h
    
    echo "uploading logs"
    gcloud storage cp /cromwell_root/stdout "$log_output_path/positioning.out"
    gcloud storage cp /cromwell_root/stderr "$log_output_path/positioning.err"
    cat stdout stderr > positioning.log
    gcloud storage cp positioning.usage "$log_output_path/positioning.usage"
    echo "<< completed positioning >>"
  >>>
 
  output {
    File coords = "coords.csv"
    File summary = "summary.pdf"
  }
  
  runtime {
    docker: docker
    memory: "~{mem_GiB} GB"
    disks: "local-disk ~{disk_GiB} SSD"
    cpu: 16
    preemptible: 0
  }

}
