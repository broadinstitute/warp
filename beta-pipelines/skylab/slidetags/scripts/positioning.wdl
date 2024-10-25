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
    set -euo pipefail
    set -x
    echo "<< starting spatial-count >>"

    gcloud config set storage/process_count 16
    gcloud config set storage/thread_count  2

    # Download the scripts
    wget https://raw.githubusercontent.com/aawdeh/Macosko-Pipelines/refs/heads/aa-pos/positioning/load_matrix.R
    wget https://raw.githubusercontent.com/aawdeh/Macosko-Pipelines/refs/heads/aa-pos/positioning/positioning.R
    wget https://raw.githubusercontent.com/aawdeh/Macosko-Pipelines/refs/heads/aa-pos/positioning/run-positioning.R

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
    ls output/* 

    if [[ -f output/seurat.qs ]] ; then
        echo "true" > DONE
    else
        echo ; echo "ERROR: CANNOT FIND: seurat.qs"
    fi

    echo; echo "Writing logs:"
    echo; echo "RNA size:"; du -sh RNA
    echo; echo "SB size:"; du -sh SB
    echo; echo "output size:"; du -sh output
    echo; echo "FREE SPACE:"; df -h
    
    echo "tar files/logs"
    cat stdout stderr > positioning.log
    tar -zcvf output.tar.gz output
    echo "<< completed positioning >>"
  >>>
 
  output {
    File output_file = "output.tar.gz"
    File positioning_log = "positioning.log"
  }
  
  runtime {
    docker: docker
    memory: "~{mem_GiB} GB"
    disks: "local-disk ~{disk_GiB} SSD"
    cpu: 16
  }

}
