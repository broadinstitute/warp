version 1.0
task PairedTagDemultiplex {
    input {
        File read1_fastq
        File read3_fastq
        File barcodes_fastq
        String input_id
        Boolean preindex
        File whitelist
        String docker_path

        Int cpu = 1
        Int disk_size = ceil(2 * (size(read1_fastq, "GiB") + size(read3_fastq, "GiB") + size(barcodes_fastq, "GiB") )) + 400
        Int preemptible = 3
        Int mem_size = 8
    }
    String r1_base = basename(read1_fastq, ".fastq.gz")
    String r2_base = basename(barcodes_fastq, ".fastq.gz")
    String r3_base = basename(read3_fastq, ".fastq.gz")
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
        docker_path: "(optional) the docker image containing the runtime environment for this task"
        mem_size: "(optional) the amount of memory (MiB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
        disk_size: "(optional) the amount of disk space (GiB) to provision for this task"
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"        
    }
    command <<<
        set -e

        # echo the basenames
        echo "r1_base is: ~{r1_base}"
        echo "r2_base is: ~{r2_base}"
        echo "r3_base is: ~{r3_base}"

        ## Need to gunzip the r1_fastq
        pass="true"
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
        echo "performing read2 length, trimming, and orientation checks" 
        if [[ $COUNT == 27 && ~{preindex} == "false" ]]
          then
          echo "Preindex is false and length is 27 bp"
          echo "Trimming last 3 bp with UPStools"
          upstools trimfq ~{input_id}_R2.fq.gz 1 24
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
            echo "Correct barcode orientation"
          else
            pass="false"
            echo "Incorrect barcode orientation"
          fi
          echo renaming files
          mv "~{input_id}_R2_trim.fq.gz" "~{r2_base}.fq.gz"
          mv "~{input_id}_R1.fq.gz" "~{r1_base}.fq.gz"
          mv "~{input_id}_R3.fq.gz" "~{r3_base}.fq.gz"

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
            echo "Correct barcode orientation"
          else
            pass="false"
            echo "Incorrect barcode orientation"
          fi
          # rename files to original name
          mv "~{input_id}_R2_prefix.fq.gz" "~{r2_base}.fq.gz"
          mv "~{input_id}_R1_prefix.fq.gz" "~{r1_base}.fq.gz"
          mv "~{input_id}_R3_prefix.fq.gz" "~{r3_base}.fq.gz"
        elif [[ $COUNT == 24 && ~{preindex} == "false" ]]
          then
          echo "FASTQ has correct index length, no modification necessary"

          ls -lh

          mv "~{input_id}_R2.fq.gz" "~{r2_base}.fq.gz"
          mv "~{input_id}_R1.fq.gz" "~{r1_base}.fq.gz"
          mv "~{input_id}_R3.fq.gz" "~{r3_base}.fq.gz"
        elif [[ $COUNT == 24 && ~{preindex} == "true" ]]
          then
          pass="false"
          echo "FASTQ does not have correct length for preindexing option"        
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
        docker: docker_path
        cpu: cpu
        memory: "${mem_size} GiB"
        disks: "local-disk ${disk_size} HDD"
        preemptible: preemptible        
    }

    output {
        File fastq1 = "~{r1_base}.fq.gz"
        File barcodes = "~{r2_base}.fq.gz"
        File fastq3 = "~{r3_base}.fq.gz"
    }
}

task AddBBTag {
    input {
        File bam
        String input_id
        String docker_path

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
        docker_path: "The docker image path containing the runtime environment for this task"
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
        docker: docker_path
        cpu: cpu
        memory: "${mem_size} GiB"
        disks: "local-disk ${disk_size} HDD"
        preemptible: preemptible
    }

    output {
        File bb_bam = "~{input_id}.bam.BB.bam"
    }
}

task ParseBarcodes {
    input {
        File atac_h5ad
        File atac_fragment
        Int nthreads = 1
        String cpuPlatform = "Intel Cascade Lake"
        String docker_path
        Int disk =  ceil((size(atac_h5ad, "GiB") + size(atac_fragment, "GiB")) * 8) + 100
        Int machine_mem_mb = ceil((size(atac_h5ad, "MiB") + size(atac_fragment, "MiB")) * 8) + 10000
    }

    String atac_base_name = basename(atac_h5ad, ".h5ad")
    String atac_fragment_base = basename(atac_fragment, ".sorted.tsv.gz")


  parameter_meta {
      atac_h5ad: "The resulting h5ad from the ATAC workflow."
      atac_fragment: "The resulting fragment TSV from the ATAC workflow."
  }

  command <<<
    set -e pipefail

    # decompress the bgzipped atac file
    echo "Moving fragment tsv for decompression" 
    mv ~{atac_fragment} ~{atac_fragment_base}.sorted.tsv.gz
    echo "Decompressing fragment file"
    bgzip -d "~{atac_fragment_base}.sorted.tsv.gz"
    echo "Done decompressing"

    python3 <<CODE

    # set parameters
    atac_h5ad = "~{atac_h5ad}"
    atac_fragment = "~{atac_fragment_base}.sorted.tsv"

    # import anndata to manipulate h5ad files
    import anndata as ad
    import pandas as pd
    import snapatac2 as snap
    print("Reading ATAC h5ad:")
    atac_data = ad.read_h5ad("~{atac_h5ad}")
    print("Reading ATAC fragment file:")
    test_fragment = pd.read_csv(atac_fragment, sep="\t", names=['chr','start', 'stop', 'barcode','n_reads'])


    # Separate out CB and preindex in the h5ad and identify sample barcodes assigned to more than one cell barcode
    print("Setting preindex and CB columns in h5ad")
    df_h5ad = atac_data.obs
    df_h5ad["preindex"] = df_h5ad.index.str[:3]
    df_h5ad["CB"] = df_h5ad.index.str[3:]
    df_h5ad["duplicates"] = 0
    preindex_counts = df_h5ad.groupby('CB')['preindex'].nunique()
    df_h5ad.loc[df_h5ad['CB'].isin(preindex_counts[preindex_counts > 1].index), 'duplicates'] = 1
      
    # Separate out CB and preindex in the fragment file
    print("Setting preindex and CB columns in fragment file")
    test_fragment["preindex"] = test_fragment["barcode"].str[:3]
    test_fragment["CB"] = test_fragment["barcode"].str[3:]
      
    # Create a new column 'duplicates' initialized with 0
    test_fragment['duplicates'] = 0
      
    # Group by 'CB' and count the number of unique 'preindex' values for each group
    preindex_counts = test_fragment.groupby('CB')['preindex'].nunique()
      
    # Update the 'duplicates' column for rows with more than one unique 'preindex' for a 'CB'
    test_fragment.loc[test_fragment['CB'].isin(preindex_counts[preindex_counts > 1].index), 'duplicates'] = 1
      
    # Idenitfy the barcodes in the whitelist that match barcodes in datasets
    atac_data.write_h5ad("~{atac_base_name}.h5ad")
    test_fragment.to_csv("~{atac_fragment_base}.tsv", sep='\t', index=False, header = False)
    CODE
    
    # sorting the file
    echo "Sorting file"
    sort -k1,1V -k2,2n "~{atac_fragment_base}.tsv" > "~{atac_fragment_base}.sorted.tsv"
    echo "Starting bgzip"
    bgzip "~{atac_fragment_base}.sorted.tsv"
    echo "Starting tabix"
    tabix -s 1 -b 2 -e 3 -C "~{atac_fragment_base}.sorted.tsv.gz"

  >>>

  runtime {
      docker: docker_path
      disks: "local-disk ~{disk} HDD"
      memory: "${machine_mem_mb} MiB"
      cpu: nthreads
  }

  output {
      File atac_h5ad_file = "~{atac_base_name}.h5ad"
      File atac_fragment_tsv = "~{atac_fragment_base}.sorted.tsv.gz"
      File atac_fragment_tsv_tbi = "~{atac_fragment_base}.sorted.tsv.gz.csi"
  }
}

task MaskPeakCallingMetrics {
  input {
    File library_metrics 
    String docker_path
  }
  command <<<
    set -e

    # Array of strings to filter out
    filter_strings=(
      "tss_enrichment_score"
      "fraction_of_high-quality_fragments_overlapping_tss"
      "fraction_of_genome_in_peaks"
      "fraction_of_high-quality_fragments_overlapping_peaks"
      "number_of_peaks"
      "fraction_of_transposition_events_in_peaks_in_cells"
    )

    # Copy input file to working directory
    cp ~{library_metrics} library_metrics_input.txt

    # Remove lines starting with any of the filter strings
    for string in "${filter_strings[@]}"; do
      grep -v "^${string}" library_metrics_input.txt > temp_file && mv temp_file library_metrics_input.txt
    done

    # Output the filtered file
    mv library_metrics_input.txt library_metrics_filtered.txt

  >>>

  runtime {
    docker: docker_path
    disks: "local-disk 32 HDD"
    memory: "8000 MiB"
  }
  output {
    File library_metrics_file = "library_metrics_filtered.txt"
  }
}