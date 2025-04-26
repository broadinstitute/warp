version 1.0

task generate_positioning {
  input {
    Array[String] rna_paths
    String sb_path
    String input_id
    Int mem_GiB  = 128
    Int disk_GiB = 128
    Int nthreads = 16
    String docker
  }
  command <<<
    set -euo pipefail
    set -x
    echo "<< starting spatial-count >>"
    
    Rscript -e "install.packages(c('optparse', 'BiocManager'), repos='https://cloud.r-project.org'); BiocManager::install('IRanges')"

    gcloud config set storage/process_count 16 # is this set by user?
    gcloud config set storage/thread_count  2 # is this set by user?

    # Download the scripts -- these need to be changed -- also need to add to docker
<<<<<<< HEAD
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/run-positioning.R
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/positioning.R
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/helpers.R   
  
=======
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/d89176cf21e072fe8b5aad3a1454ad194fca7c9a/slide-tags/run-positioning.R
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/d89176cf21e072fe8b5aad3a1454ad194fca7c9a/slide-tags/positioning.R
    wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/d89176cf21e072fe8b5aad3a1454ad194fca7c9a/slide-tags/helpers.R

>>>>>>> ff3e3e0fbebf9d7c65ee0f7242007ec4cb7ae89e
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
    gcloud storage cp ~{sb_path} .
    baseSB=`basename ~{sb_path}`

    # Run the script
    echo ; echo "Running run-positioning.R"
    Rscript run-positioning.R RNA $baseSB output

    # Upload the results
    ls output/* 

    if [[ -f output/seurat.qs ]] ; then
        echo "true" > DONE
    else
        echo ; echo "ERROR: CANNOT FIND: seurat.qs"
    fi

    echo; echo "Writing logs:"
    echo; echo "RNA size:"; du -sh RNA
    echo; echo "SB size:"; du -sh $baseSB
    echo; echo "output size:"; du -sh output
    echo; echo "FREE SPACE:"; df -h
    
    echo "tar files/logs"
    cat stdout stderr > positioning.log
    
    # Rename and move files
    mv output/* .
    mv summary.pdf ~{input_id}_summary.pdf
    mv seurat.qs ~{input_id}_seurat.qs
    mv coords.csv ~{input_id}_coords.csv
<<<<<<< HEAD

    tar -zcvf output.tar.gz output/matrix.csv.gz output/cb_whitelist.txt output/spatial_metadata.json
    mv output.tar.gz ~{input_id}_output.tar.gz
=======
    mv coords2.csv ~{input_id}_coords2.csv
    
    ls 
    tar -zcvf output.tar.gz matrix.csv.gz cb_whitelist.txt spatial_metadata.json
    mv output.tar.gz ~{input_id}_intermediates.tar.gz
>>>>>>> ff3e3e0fbebf9d7c65ee0f7242007ec4cb7ae89e
    mv positioning.log ~{input_id}_positioning.log
    echo "<< completed positioning >>"
  >>>
 
  output {
    File seurat_qs = "~{input_id}_seurat.qs"
    File coords_csv = "~{input_id}_coords.csv"
<<<<<<< HEAD
    File summary_pdf = "~{input_id}_summary.pdf"
    File output_file = "~{input_id}_output.tar.gz"
=======
    File coords2_csv = "~{input_id}_coords2.csv"
    File summary_pdf = "~{input_id}_summary.pdf"
    File intermediates_file = "~{input_id}_intermediates.tar.gz"
>>>>>>> ff3e3e0fbebf9d7c65ee0f7242007ec4cb7ae89e
    File positioning_log = "~{input_id}_positioning.log"
  }
  
  runtime {
    docker: docker
    memory: "~{mem_GiB} GB"
    disks: "local-disk ~{disk_GiB} SSD"
    cpu: nthreads
  }

}
