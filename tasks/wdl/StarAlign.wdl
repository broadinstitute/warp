version 1.0

task StarAlignBamSingleEnd {
  input {
    File bam_input
    File tar_star_reference

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/star:1.0.0-2.7.9a-1658781884"
    Int machine_mem_mb = ceil((size(tar_star_reference, "Gi")) + 6) * 1100
    Int cpu = 16
    # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
    Int disk = ceil((size(tar_star_reference, "Gi") * 2.5) + (size(bam_input, "Gi") * 2.5))
    # by default request non preemptible machine to make sure the slow star alignment step completes
    Int preemptible = 0
  }

  meta {
    description: "Aligns reads in bam_input to the reference genome in tar_star_reference"
  }

  parameter_meta {
    bam_input: "unaligned bam file containing genomic sequence, tagged with barcode information"
    tar_star_reference: "star reference tarball built against the species that the bam_input is derived from"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    # prepare reference
    mkdir genome_reference
    tar -xf "${tar_star_reference}" -C genome_reference --strip-components 1
    rm "${tar_star_reference}"

    STAR \
      --runMode alignReads \
      --runThreadN ${cpu} \
      --genomeDir genome_reference \
      --readFilesIn "${bam_input}" \
      --outSAMtype BAM Unsorted \
      --outSAMmultNmax -1 \
      --outSAMattributes All \
      --outSAMunmapped Within \
      --readFilesType SAM SE \
      --readFilesCommand samtools view -h \
      --runRNGseed 777
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} SSD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File bam_output = "Aligned.out.bam"
    File alignment_log = "Log.final.out"
  }
}

task StarAlignFastqPairedEnd {
  input {
    File fastq1
    File fastq2
    File tar_star_reference

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/star:1.0.0-2.7.9a-1658781884"
    Int machine_mem_mb = ceil((size(tar_star_reference, "Gi")) + 6) * 1100
    Int cpu = 16
    # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
    Int disk = ceil((size(tar_star_reference, "Gi") * 2.5) + (size(fastq1, "Gi") * 5.0))
    # by default request non preemptible machine to make sure the slow star alignment step completes
    Int preemptible = 3
  }

  meta {
    description: "Aligns reads in fastq1 and fastq2 to the reference genome in tar_star_reference"
  }

  parameter_meta {
    fastq1: "trimmed R1 FASTQ file containing genomic sequence"
    fastq2: "trimmed R2 FASTQ file containing genomic sequence"
    tar_star_reference: "star reference tarball built against the species that the bam_input is derived from"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    # prepare reference
    mkdir genome_reference
    tar -xf "${tar_star_reference}" -C genome_reference --strip-components 1
    rm "${tar_star_reference}"

    STAR \
    --genomeDir genome_reference \
    --runThreadN ${cpu} \
    --readFilesIn ~{fastq1} ~{fastq2} \
    --readFilesCommand "gunzip -c" \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx \
    --runRNGseed 777 \
    --limitBAMsortRAM 10000000000 \
    --quantMode GeneCounts
  }
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} SSD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File output_bam = "Aligned.sortedByCoord.out.bam"
  }

}

task StarAlignFastqMultisample {
  input {
    Array[File] fastq1_input_files
    Array[File] fastq2_input_files
    Array[String] input_ids
    File tar_star_reference

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/star:1.0.0-2.7.9a-1658781884"
    Int machine_mem_mb = ceil((size(tar_star_reference, "Gi")) + 6) * 1100
    Int cpu = 16
    # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
    Int disk = ceil((size(tar_star_reference, "Gi") * 2.5) + (size(fastq2_input_files, "Gi") * 2.0))
    # by default request non preemptible machine to make sure the slow star alignment step completes
    Int preemptible = 3
  }

  meta {
    description: "Aligns reads in fastq1 and fastq2 to the reference genome in tar_star_reference"
  }

  parameter_meta {
    fastq1_input_files: "Array of trimmed R1 fastq files containing genomic sequence."
    fastq2_input_files: "Array of trimmed R2 fastq files containing genomic sequence."
    input_ids: "Array of input ids"
    tar_star_reference: "star reference tarball built against the species that the bam_input is derived from"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command <<<
    set -e
    set -exo pipefail

    # prepare reference
    mkdir genome_reference
    tar -xf "~{tar_star_reference}" -C genome_reference --strip-components 1
    rm "~{tar_star_reference}"

    declare -a fastq1_files=(~{sep=' ' fastq1_input_files})
    declare -a fastq2_files=(~{sep=' ' fastq2_input_files})
    declare -a output_prefix=(~{sep=' ' input_ids})
    STAR --genomeLoad LoadAndExit --genomeDir genome_reference
    for (( i=0; i<${#output_prefix[@]}; ++i));
      do
        STAR \
          --genomeDir genome_reference \
          --runThreadN ~{cpu} \
          --readFilesIn ${fastq1_files[$i]} ${fastq2_files[$i]} \
          --readFilesCommand "gunzip -c" \
          --outSAMtype BAM SortedByCoordinate \
          --outReadsUnmapped Fastx \
          --runRNGseed 777 \
          --limitBAMsortRAM 10000000000 \
          --quantMode GeneCounts \
          --genomeLoad LoadAndKeep

        mv "Aligned.sortedByCoord.out.bam"   "${output_prefix[$i]}_Aligned.bam"
      done;
    STAR --genomeLoad Remove --genomeDir genome_reference
  >>>

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} SSD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }

  output {
    Array[File] output_bam = glob("*_Aligned.bam")
  }

}


task STARsoloFastq {
  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    File tar_star_reference
    File white_list
    Int chemistry
    String star_strand_mode
    String counting_mode # when counting_mode = sn_rna, runs Gene and GeneFullEx50pAS in single alignments
    String input_id
    String output_bam_basename
    Boolean? count_exons
    String? soloMultiMappers
    String soloCBmatchWLtype = "1MM_multi" #"1MM_multi_Nbase_pseudocounts"
    Int expected_cells = 3000
    String reference_path = tar_star_reference

    # runtime values
    String samtools_star_docker_path
    String cpu_platform = "Intel Ice Lake"
    Int input_size = ceil(size(r1_fastq, "GiB") + size(r2_fastq, "GiB"))
    Int cpu = 16
    Int disk = 5000
    Int limitBAMsortRAM = 30
    Int mem_size = if ceil(size(r1_fastq, "GiB") + size(r2_fastq, "GiB")) <= 100 then 64 else 128

    # by default request non preemptible machine to make sure the slow star alignment step completes
    Int preemptible = 1
  }

  Int outBAMsortingBinsN = (((ceil(size(r1_fastq, "GiB") + size(r2_fastq, "GiB")) + 50) / 100) * 100) + 100

  meta {
    description: "Aligns reads in bam_input to the reference genome in tar_star_reference" 
  }

  parameter_meta {
    r1_fastq: "input FASTQ file array"
    r2_fastq: "array of forward read FASTQ files"
    tar_star_reference: "star reference tarball built against the species that the bam_input is derived from"
    star_strand_mode: "STAR mode for handling stranded reads. Options are 'Forward', 'Reverse, or 'Unstranded'"
    samtools_star_docker_path: "(optional) the docker image containing the runtime environment for this task"
    mem_size: "the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    limitBAMsortRAM: "(optional) Specifies the maximum amount of RAM (in GiB) allocated for sorting BAM files in STAR. Default is 30."
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command <<<
    set -euo pipefail
    set -x

    ulimit -n 10000

    UMILen=10
    CBLen=16
    if [ "~{chemistry}" == 2 ]
    then
        ## V2
        UMILen=10
        CBLen=16
    elif [ "~{chemistry}" == 3 ]
    then
        ## V3
        UMILen=12
        CBLen=16
    else
        echo Error: unknown chemistry value: "$chemistry". Should be one of "tenX_v2" or "texX_v3".
        exit 1;
    fi

    # Check that the star strand mode matches STARsolo aligner options
    if [[ "~{star_strand_mode}" == "Forward" ]] || [[ "~{star_strand_mode}" == "Reverse" ]] || [[ "~{star_strand_mode}" == "Unstranded" ]]
    then
        ## single cell or whole cell
        echo STAR mode is assigned
    else
        echo Error: unknown STAR strand mode: "~{star_strand_mode}". Should be Forward, Reverse, or Unstranded.
        exit 1;
    fi

    # prepare reference
    mkdir genome_reference
    tar -xf "~{tar_star_reference}" -C genome_reference --strip-components 1
    rm "~{tar_star_reference}"

    COUNTING_MODE=""
    if [[ "~{counting_mode}" == "sc_rna" ]]
    then
        # single cell or whole cell
        COUNTING_MODE="Gene"
        echo "Running in ~{counting_mode} mode. The Star parameter --soloFeatures will be set to $COUNTING_MODE"
    elif [[ "~{counting_mode}" == "sn_rna" ]]
    then
        # single nuclei
        if [[ ~{count_exons} == false ]]
        then
            COUNTING_MODE="GeneFull_Ex50pAS"
            echo "Running in ~{counting_mode} mode. Count_exons is false and the Star parameter --soloFeatures will be set to $COUNTING_MODE"
        else
            COUNTING_MODE="GeneFull_Ex50pAS Gene"
            echo "Running in ~{counting_mode} mode. Count_exons is true and the Star parameter --soloFeatures will be set to $COUNTING_MODE"     
        fi
    else
        echo Error: unknown counting mode: "$counting_mode". Should be either sn_rna or sc_rna.
        exit 1;
    fi

    # convert limitBAMsortRAM from GB to bytes 
    RAM_limit_bytes=$((1073741824 * ~{limitBAMsortRAM})) 
    echo $RAM_limit_bytes, ~{limitBAMsortRAM}
    
    # run star
    STAR \
        --soloType Droplet \
        --soloStrand ~{star_strand_mode} \
        --runThreadN ~{cpu} \
        --genomeDir genome_reference \
        --readFilesIn "~{sep=',' r2_fastq}" "~{sep=',' r1_fastq}" \
        --readFilesCommand "gunzip -c" \
        --soloCBwhitelist ~{white_list} \
        --soloUMIlen $UMILen --soloCBlen $CBLen \
        --soloFeatures $COUNTING_MODE \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30  \
        --soloCBmatchWLtype ~{soloCBmatchWLtype} \
        --soloUMIdedup 1MM_CR \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes UB UR UY CR CB CY NH GX GN sF cN \
        --soloBarcodeReadLength 0 \
        --limitBAMsortRAM $RAM_limit_bytes \
        --outBAMsortingBinsN ~{outBAMsortingBinsN} \
        --soloCellReadStats Standard \
        ~{"--soloMultiMappers " + soloMultiMappers} \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloCellFilter EmptyDrops_CR

    # validate the bam with samtools quickcheck
    samtools quickcheck -v Aligned.sortedByCoord.out.bam
    # reheader the BAM
    samtools view -H Aligned.sortedByCoord.out.bam > header.txt
    echo -e "@CO\tReference genome used: ~{reference_path}" >> header.txt
    samtools reheader header.txt Aligned.sortedByCoord.out.bam > Aligned.sortedByCoord.out.reheader.bam

    echo "UMI LEN " $UMILen
    touch barcodes_sn_rna.tsv features_sn_rna.tsv matrix_sn_rna.mtx CellReads_sn_rna.stats Features_sn_rna.stats Summary_sn_rna.csv UMIperCellSorted_sn_rna.txt
      
    ###########################################################################
    # SAVE OUTPUT FILES
    ###########################################################################
    # Function to move .mtx files to /cromwell_root/
    move_mtx_files() {
      local directory=$1
      echo "Processing $directory"
      find "${directory}/raw" -maxdepth 1 -type f -name "*.mtx" -print0 | xargs -0 -I{} sh -c 'echo Moving {}; mv {} /cromwell_root/'
    }

    # Function to move and rename common files
    move_common_files() {
      local src_dir=$1
      local suffix=$2

      declare -A files=(
            ["barcodes.tsv"]="barcodes.tsv"
            ["features.tsv"]="features.tsv"
            ["CellReads.stats"]="CellReads.stats"
            ["Features.stats"]="Features.stats"
            ["Summary.csv"]="Summary.csv"
            ["UMIperCellSorted.txt"]="UMIperCellSorted.txt"
      )

      for file in "${!files[@]}"; do
          file_path="${files[$file]}"
          name=$(basename "$file_path")
          base="${name%.*}"
          extension="${name##*.}"
          new_name="${base}${suffix}.${extension}"
          echo $new_name
          if [[ -f "$src_dir/raw/$file" ]]; then
                echo "Renaming $src_dir/raw/$file → $new_name"
                mv "$src_dir/raw/$file" "$new_name"
          elif [[ -f "$src_dir/$file" ]]; then
                echo "Renaming $src_dir/$file → $new_name"
                mv "$src_dir/$file" "$new_name"
          else
                echo "Warning: Missing file in $src_dir or $src_dir/raw: $file"
          fi
      done
    }

    if [[ "~{counting_mode}" == "sc_rna" ]]
    then
      SoloDirectory="Solo.out/Gene"
      echo "SoloDirectory is $SoloDirectory"
      move_mtx_files "$SoloDirectory"
      move_common_files "$SoloDirectory" ""
    elif [[ "~{counting_mode}" == "sn_rna" ]]
    then
      SoloDirectory="Solo.out/GeneFull_Ex50pAS"
      move_mtx_files "$SoloDirectory"     
      if [[ "~{count_exons}" == "true" ]]; then
        # Additional processing for sn_rna with exon counting
        SoloDirectory2="Solo.out/Gene"
        find "$SoloDirectory2/raw" -maxdepth 1 -type f -name "*.mtx" -print0 | xargs -0 -I{} sh -c 'new_name="$(basename {} .mtx)_sn_rna.mtx"; echo Renaming {}; mv {} "/cromwell_root/$new_name"'
        move_common_files "$SoloDirectory2" "_sn_rna"  # Add snRNA for renaming
      fi
      move_common_files "$SoloDirectory" ""  # Standard snRNA renaming
    else
      echo Error: unknown counting mode: "$counting_mode". Should be either sn_rna or sc_rna.
    fi

    # filtered outputs in Solo.out/GeneFull_Ex50pAS/filtered: barcodes.tsv features.tsv matrix.mtx
    ls ${SoloDirectory}/filtered
    echo "Tarring up filtered matrix files"
    tar -cvf ~{input_id}_filtered_mtx_files.tar ${SoloDirectory}/filtered/barcodes.tsv ${SoloDirectory}/filtered/features.tsv ${SoloDirectory}/filtered/matrix.mtx
    echo "Done processing"

    # List the final directory contents
    echo "Final directory listing:"
    ls -l
    mv Aligned.sortedByCoord.out.reheader.bam ~{output_bam_basename}.bam

    ###########################################################################
    # FROM MERGE STAR OUTPUT TASK
    ###########################################################################
    # Function to process a matrix (regular or snRNA)
    process_matrix() {
        local MATRIX_NAME=$1  # matrix or matrix_sn_rna
        local BARCODE_FILE=$2
        local FEATURE_FILE=$3
        local MATRIX_FILE=$4
        local OUTPUT_DIR=$5

        echo "Processing $MATRIX_NAME data..."

        # Create and copy matrix files
        mkdir -p ./$MATRIX_NAME
        cp $MATRIX_FILE ./$MATRIX_NAME/matrix.mtx
        cp $BARCODE_FILE ./$MATRIX_NAME/barcodes.tsv
        cp $FEATURE_FILE ./$MATRIX_NAME/features.tsv

        # Compress matrix files
        tar -zcvf ~{input_id}_${MATRIX_NAME}.mtx_files.tar -C ./$MATRIX_NAME .

        # List files
        echo "Listing files after processing $MATRIX_NAME:"
        ls

        # If text files are present, create a tar archive with them and run python script to combine shard metrics
        python3 /scripts/scripts/combine_shard_metrics.py \
          Summary.csv Features.stats CellReads.stats ~{counting_mode} ~{input_id} ${SoloDirectory}/filtered/barcodes.tsv ${SoloDirectory}/filtered/matrix.mtx ~{expected_cells}

        echo "tarring STAR txt files"
        tar -zcvf ~{input_id}.star_metrics.tar *.txt
       
        # Create the compressed raw count matrix
        python3 /scripts/scripts/create-merged-npz-output.py \
            --barcodes $BARCODE_FILE --features $FEATURE_FILE --matrix $MATRIX_FILE --input_id ~{input_id}
     
      }

    # Process main matrix
    process_matrix "matrix" "barcodes.tsv" "features.tsv" "matrix.mtx" "./output"

    # Process snRNA matrix only if files exist
    if [ -s "barcodes_sn_rna.tsv" ]; then
        process_matrix "matrix_sn_rna" "barcodes_sn_rna.tsv" "features_sn_rna.tsv" "matrix_sn_rna.mtx" "./outputsnrna"
    fi

    ls -lR

    mv ${SoloDirectory}/filtered/barcodes.tsv filtered_barcodes.tsv
    cat filtered_barcodes.tsv
  >>>

  runtime {
    docker: samtools_star_docker_path
    memory: "~{mem_size} GiB"
    disks: "local-disk ~{disk} SSD"
    disk: disk + " GB" # TES
    cpu: cpu
    cpuPlatform: cpu_platform
    preemptible: preemptible
  }

  output {
    File bam_output = "~{output_bam_basename}.bam"
    File alignment_log = "Log.final.out"
    File general_log = "Log.out"
    File barcodes = "barcodes.tsv"
    File features = "features.tsv"
    File matrix = "matrix.mtx"
    File barcodes_sn_rna = "barcodes_sn_rna.tsv"
    File features_sn_rna = "features_sn_rna.tsv"
    File matrix_sn_rna = "matrix_sn_rna.mtx"
    File cell_reads = "CellReads.stats"
    File align_features = "Features.stats"
    File summary = "Summary.csv"
    File umipercell = "UMIperCellSorted.txt"
    File cell_reads_sn_rna = "CellReads_sn_rna.stats"
    File align_features_sn_rna = "Features_sn_rna.stats"
    File summary_sn_rna = "Summary_sn_rna.csv"
    File umipercell_sn_rna = "UMIperCellSorted_sn_rna.txt"
    File? multimappers_EM_matrix = "UniqueAndMult-EM.mtx"
    File? multimappers_Uniform_matrix = "UniqueAndMult-Uniform.mtx"
    File? multimappers_Rescue_matrix = "UniqueAndMult-Rescue.mtx"
    File? multimappers_PropUnique_matrix = "UniqueAndMult-PropUnique.mtx"
    # Output files for previous merging STAR output step
    File row_index = "~{input_id}_sparse_counts_row_index.npy"
    File col_index = "~{input_id}_sparse_counts_col_index.npy"
    File sparse_counts = "~{input_id}_sparse_counts.npz"
    File? library_metrics="~{input_id}_library_metrics.csv"
    File? mtx_files ="~{input_id}.mtx_files.tar"
    File? filtered_mtx_files = "~{input_id}_filtered_mtx_files.tar"
    File? cell_reads_out = "~{input_id}.star_metrics.tar"
    File outputbarcodes = "filtered_barcodes.tsv"
  }
}

task MergeStarOutput {

  input {
    Array[File] barcodes
    Array[File] features
    Array[File] matrix
    Array[File]? cell_reads
    Array[File]? summary
    Array[File]? align_features
    Array[File]? umipercell
    String? counting_mode
    
    String input_id
    # additional library aliquot id
    String gex_nhash_id = ""
    Int expected_cells = 3000
    File barcodes_single = barcodes[0]
    File features_single = features[0]

    #runtime values
    String star_merge_docker_path

    Int machine_mem_gb = 20
    Int cpu = 1
    Int disk = ceil(size(matrix, "Gi") * 2) + 10
    Int preemptible = 3
  }
  meta {
    description: "Create three files as .npy,  .npy and .npz for numpy array and scipy csr matrix for the barcodes, gene names and the count matrix by merging multiple  STARSolo count matrices (mtx format)."
  }

  parameter_meta {
    star_merge_docker_path: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_gb: "(optional) the amount of memory (GiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command <<<
    set -euo pipefail
    set -x 

    declare -a barcodes_files=(~{sep=' ' barcodes})
    declare -a features_files=(~{sep=' ' features})
    declare -a matrix_files=(~{sep=' ' matrix})
    declare -a cell_reads_files=(~{sep=' ' cell_reads})
    declare -a summary_files=(~{sep=' ' summary})
    declare -a align_features_files=(~{sep=' ' align_features})
    declare -a umipercell_files=(~{sep=' ' umipercell})

    # create the  compressed raw count matrix with the counts, gene names and the barcodes
    python3 /scripts/scripts/combined_mtx.py \
    ${matrix_files[@]} \
    ~{input_id}.uniform.mtx

    mkdir matrix
    #Using cp because mv isn't moving
    pwd
    ls -lR
    cp ~{input_id}.uniform.mtx ./matrix/matrix.mtx
    cp ~{barcodes_single} ./matrix/barcodes.tsv
    cp ~{features_single} ./matrix/features.tsv

    tar -zcvf ~{input_id}.mtx_files.tar ./matrix/*


    # Running star for combined cell matrix
    # outputs will be called outputbarcodes.tsv. outputmatrix.mtx, and outputfeatures.tsv
    STAR --runMode soloCellFiltering ./matrix ./output --soloCellFilter EmptyDrops_CR
    
    #list files
    echo "listing files"
    ls
    # if theres a file in cell_reads_files --  check if non empty
    if [ -n "${cell_reads_files[*]}" ]; then
      # Destination file for cell reads
      dest="~{input_id}_cell_reads.txt"
    
      # first create the header from the first file in the list, and add a column header for the shard id
      head -n 1 "${cell_reads_files[0]}" | awk '{print $0 "\tshard_number"}' > "$dest"
    
      # Loop through the array and add the second row with shard number to a temp file notinpasslist.txt
      for index in "${!cell_reads_files[@]}"; do
        secondLine=$(sed -n '2p' "${cell_reads_files[$index]}")
        echo -e "$secondLine\t$index" >> "notinpasslist.txt"
      done

      # add notinpasslist.txt to the destination file and delete the notinpasslist.txt
      cat "notinpasslist.txt" >> "$dest"
      rm notinpasslist.txt

      # now add the shard id to the matrix in a temporary matrix file, and skip the first two lines
      counter=0
      for cell_read in "${cell_reads_files[@]}"; do
        if [ -f "$cell_read" ]; then
          awk -v var="$counter" 'NR>2 {print $0 "\t" var}' "$cell_read" >> "matrix.txt" 
          let counter=counter+1
        fi
      done

      # add the matrix to the destination file, then delete the matrix file
      cat "matrix.txt" >> "$dest"
      rm "matrix.txt"
    fi

    counter=0
    for summary in "${summary_files[@]}"; do
      if [ -f "$summary" ]; then
        awk -v var=",$counter" '{print $0 var}' "$summary" >> "~{input_id}_summary.txt"
        let counter=counter+1
      fi
    done
    
    counter=0
    for align_feature in "${align_features_files[@]}"; do
      if [ -f "$align_feature" ]; then
        awk -v var="$counter" '{print $0 " " var}' "$align_feature" >> "~{input_id}_align_features.txt"
        let counter=counter+1
      fi
    done

    # note that the counter might not correspond to the shard number, it is just the order of files in bash (e.g. 10 before 2)
    counter=0
    for umipercell in "${umipercell_files[@]}"; do
      if [ -f "$umipercell" ]; then
        awk -v var="$counter" '{print $0, var}' "$umipercell" >> "~{input_id}_umipercell.txt"
        let counter=counter+1
      fi
    done
    
    # If text files are present, create a tar archive with them and run python script to combine shard metrics
    if ls *.txt 1> /dev/null 2>&1; then
      echo "listing files"
      ls
      python3 /scripts/scripts/combine_shard_metrics.py \
      ~{input_id}_summary.txt \
      ~{input_id}_align_features.txt \
      ~{input_id}_cell_reads.txt \
      ~{counting_mode} \
      ~{input_id} \
      outputbarcodes.tsv \
      outputmatrix.mtx \
      ~{expected_cells}

      echo "tarring STAR txt files"
      tar -zcvf ~{input_id}.star_metrics.tar *.txt
    else
      echo "No text files found in the folder."
    fi

   #
   # create the  compressed raw count matrix with the counts, gene names and the barcodes
    python3 /scripts/scripts/create-merged-npz-output.py \
        --barcodes ${barcodes_files[@]} \
        --features ${features_files[@]} \
        --matrix ${matrix_files[@]} \
        --input_id ~{input_id}

    # tar up filtered matrix outputbarcodes.tsv, outputfeatures.tsv, outputmatrix.mtx
    echo "Tarring up filtered matrix files"
    tar -cvf ~{input_id}_filtered_mtx_files.tar outputbarcodes.tsv outputfeatures.tsv outputmatrix.mtx
    echo "Done"
  >>>

  runtime {
    docker: star_merge_docker_path
    memory: "${machine_mem_gb} GiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File row_index = "~{input_id}_sparse_counts_row_index.npy"
    File col_index = "~{input_id}_sparse_counts_col_index.npy"
    File sparse_counts = "~{input_id}_sparse_counts.npz"
    File? cell_reads_out = "~{input_id}.star_metrics.tar"
    File? library_metrics="~{input_id}_library_metrics.csv"
    File? mtx_files ="~{input_id}.mtx_files.tar"
    File? filtered_mtx_files = "~{input_id}_filtered_mtx_files.tar"
    File? outputbarcodes = "outputbarcodes.tsv"
  }
}

task STARsoloFastqSlideSeq {
  input {
    Array[File] r1_fastq
    Array[File] r2_fastq
    File tar_star_reference
    File whitelist
    String output_bam_basename
    String read_structure
    Boolean? count_exons

    # runtime values
    String docker = "us.gcr.io/broad-gotc-prod/star:1.0.1-2.7.11a-1692706072"
    Int machine_mem_mb = 64000
    Int cpu = 8
    # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
    Int disk = ceil((size(tar_star_reference, "Gi") * 3)) + ceil(size(r1_fastq, "Gi") * 20) +  ceil(size(r2_fastq, "Gi") * 20)
    # by default request non preemptible machine to make sure the slow star alignment step completes
    Int preemptible = 3
  }

  command <<<
    set -e
    declare -a fastq1_files=(~{sep=' ' r1_fastq})
    declare -a fastq2_files=(~{sep=' ' r2_fastq})
    cut -f 1 ~{whitelist} > WhiteList.txt

    nums=$(echo ~{read_structure} | sed 's/[[:alpha:]]/ /g')
    read -a arr_num <<< $nums

    chars=$(echo ~{read_structure} | sed 's/[[:digit:]]/ /g')
    read -a arr_char <<< $chars

    UMILen=0
    CBLen=0
    for (( i=0; i<${#arr_char[@]}; ++i));
      do
        if [[ ${arr_char[$i]} == 'C' ]]
        then
          CBLen=$(( CBLen + arr_num[$i] ))
        elif [[ ${arr_char[$i]} == 'M' ]]
        then
          UMILen=$(( UMILen + arr_num[$i] ))
        fi
    done;
    UMIstart=$(( 1 + CBLen))

    # If this argument is true, we will count reads aligned to exons in addition
    COUNTING_MODE="GeneFull"
    if ~{count_exons}
    then
      COUNTING_MODE="Gene GeneFull"
    fi

    # prepare reference
    mkdir genome_reference
    tar -xf "~{tar_star_reference}" -C genome_reference --strip-components 1
    rm "~{tar_star_reference}"

    STAR \
      --soloType Droplet \
      --soloCBwhitelist WhiteList.txt \
      --soloFeatures $COUNTING_MODE \
      --runThreadN ~{cpu} \
      --genomeDir genome_reference \
      --readFilesIn $fastq2_files $fastq1_files \
      --readFilesCommand "gunzip -c" \
      --soloInputSAMattrBarcodeSeq CR UR \
      --soloInputSAMattrBarcodeQual CY UY \
      --soloCBlen $CBLen \
      --soloCBstart 1 \
      --soloUMIlen $UMILen \
      --soloUMIstart $UMIstart \
      --outSAMtype BAM SortedByCoordinate \
      --clip3pAdapterSeq AAAAAA \
      --clip3pAdapterMMp 0.1 \
      --outSAMattributes UB UR UY CR CB CY NH GX GN sF

    touch barcodes_exon.tsv
    touch features_exon.tsv
    touch matrix_exon.mtx

    mv "Solo.out/GeneFull/raw/barcodes.tsv" barcodes.tsv
    mv "Solo.out/GeneFull/raw/features.tsv" features.tsv
    mv "Solo.out/GeneFull/raw/matrix.mtx"   matrix.mtx

    if  ~{count_exons}
    then
      mv "Solo.out/Gene/raw/barcodes.tsv"     barcodes_exon.tsv
      mv "Solo.out/Gene/raw/features.tsv"     features_exon.tsv
      mv "Solo.out/Gene/raw/matrix.mtx"       matrix_exon.mtx
    fi

    mv Aligned.sortedByCoord.out.bam ~{output_bam_basename}.bam

  >>>

  runtime {
    docker: docker
    memory: "~{machine_mem_mb} MiB"
    disks: "local-disk ~{disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File bam_output = "~{output_bam_basename}.bam"
    File alignment_log = "Log.final.out"
    File general_log = "Log.out"
    File barcodes = "barcodes.tsv"
    File features = "features.tsv"
    File matrix = "matrix.mtx"
    File barcodes_sn_rna = "barcodes_exon.tsv"
    File features_sn_rna = "features_exon.tsv"
    File matrix_sn_rna = "matrix_exon.mtx"  
  }
}

task STARGenomeRefVersion {
  input {
    String tar_star_reference
    Int disk = 10
    String ubuntu_docker_path
  }

  meta {
    description: "Reads the reference file name and outputs a txt file containing genome source, build, and annotation version"
  }

  parameter_meta {
    tar_star_reference: "input STAR reference in TAR format"
  }

  command <<<
    # check genomic reference version and print to output txt file
    STRING=~{tar_star_reference}
    BASE=$(basename $STRING .tar)
    IFS=' -' read -r -a array <<< $BASE
    REFERENCE=${array[2]}
    VERSION=${array[4]}
    ANNOTATION=${array[5]}

    echo -e "$REFERENCE\n$VERSION\n$ANNOTATION" > reference_version.txt
    echo Reference is $REFERENCE
    echo Version is $VERSION
    echo Annotation is $ANNOTATION

  >>>

# Output is TXT file containing reference source, build version and annotation version
  output {
    File genomic_ref_version = "reference_version.txt"
  }

  runtime {
    docker: ubuntu_docker_path
    memory: "2 GiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu:"1"
  }
}