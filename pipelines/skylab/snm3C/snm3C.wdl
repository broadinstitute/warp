version 1.0
import "../../../tasks/broad/Utilities.wdl" as utils


workflow snm3C {

    input {
        Array[File] fastq_input_read1
        Array[File] fastq_input_read2
        File random_primer_indexes
        String plate_id
        String cloud_provider
        # mapping inputs
        File tarred_index_files
        File genome_fa
        File chromosome_sizes

        String r1_adapter = "AGATCGGAAGAGCACACGTCTGAAC"
        String r2_adapter = "AGATCGGAAGAGCGTCGTGTAGGGA"
        Int r1_left_cut = 10
        Int r1_right_cut = 10
        Int r2_left_cut = 10
        Int r2_right_cut = 10
        Int min_read_length = 30
        Int num_upstr_bases = 0
        Int num_downstr_bases = 2
        Int compress_level = 5
        Int batch_number = 6
    }
    #docker images
    String m3c_yap_hisat_docker = "m3c-yap-hisat:2.4"
    # Determine docker prefix based on cloud provider
    String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
    String acr_docker_prefix = "dsppipelinedev.azurecr.io/"
    String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix
    String cromwell_root_dir = if cloud_provider == "gcp" then "/mnt/disks/cromwell_root" else "/cromwell-executions"

    # make sure either gcp or azr is supplied as cloud_provider input
    if ((cloud_provider != "gcp") && (cloud_provider != "azure")) {
        call utils.ErrorWithMessage as ErrorMessageIncorrectInput {
        input:
            message = "cloud_provider must be supplied with either 'gcp' or 'azure'."
        }
    }

    # version of the pipeline

    String pipeline_version = "4.1.0"

    call Demultiplexing {
        input:
            fastq_input_read1 = fastq_input_read1,
            fastq_input_read2 = fastq_input_read2,
            random_primer_indexes = random_primer_indexes,
            plate_id = plate_id,
            batch_number = batch_number,
            docker = docker_prefix + m3c_yap_hisat_docker,
    }

    scatter(tar in Demultiplexing.tarred_demultiplexed_fastqs) {
        call Hisat_paired_end as Hisat_paired_end {
          input:
                tarred_demultiplexed_fastqs = tar,
                tarred_index_files = tarred_index_files,
                genome_fa = genome_fa,
                chromosome_sizes = chromosome_sizes,
                min_read_length = min_read_length,
                r1_adapter = r1_adapter,
                r2_adapter = r2_adapter,
                r1_left_cut = r1_left_cut,
                r1_right_cut = r1_right_cut,
                r2_left_cut = r2_left_cut,
                r2_right_cut = r2_right_cut,
                plate_id = plate_id,
                docker = docker_prefix + m3c_yap_hisat_docker,
                cromwell_root_dir = cromwell_root_dir,
                cloud_provider = cloud_provider,
        }

        call Hisat_single_end as Hisat_single_end {
            input:
                split_fq_tar = Hisat_paired_end.split_fq_tar,
                tarred_index_files = tarred_index_files,
                genome_fa = genome_fa,
                plate_id = plate_id,
                docker = docker_prefix + m3c_yap_hisat_docker,
                cromwell_root_dir = cromwell_root_dir,
                cloud_provider = cloud_provider
        }
        call Merge_sort_analyze as Merge_sort_analyze {
            input:
               paired_end_unique_tar = Hisat_paired_end.unique_bam_tar,
               read_overlap_tar = Hisat_single_end.remove_overlaps_output_bam_tar,
               genome_fa = genome_fa,
               num_upstr_bases = num_upstr_bases,
               num_downstr_bases = num_downstr_bases,
               compress_level = compress_level,
               chromosome_sizes = chromosome_sizes,
               plate_id = plate_id,
               docker = docker_prefix + m3c_yap_hisat_docker,
               cromwell_root_dir = cromwell_root_dir,
               cloud_provider = cloud_provider
        }
    }

    call Summary_PerCellOutput {
       input:
            name_sorted_bams = Merge_sort_analyze.name_sorted_bam,
            unique_reads_cgn_extraction_allc = Merge_sort_analyze.allc,
            unique_reads_cgn_extraction_tbi = Merge_sort_analyze.tbi,
            all_reads_3C_contacts = Merge_sort_analyze.all_reads_3C_contacts,
            unique_reads_cgn_extraction_allc_extract = Merge_sort_analyze.extract_allc_output_allc_tar,
            unique_reads_cgn_extraction_tbi_extract = Merge_sort_analyze.extract_allc_output_tbi_tar,
            plate_id = plate_id,
            docker = docker_prefix + m3c_yap_hisat_docker,
            cromwell_root_dir = cromwell_root_dir
    }

    call Summary {
        input:
            trimmed_stats = Hisat_paired_end.trim_stats_tar,
            hisat3n_stats = Hisat_paired_end.hisat3n_paired_end_stats_tar,
            r1_hisat3n_stats = Hisat_single_end.hisat3n_dna_split_reads_summary_R1_tar,
            r2_hisat3n_stats = Hisat_single_end.hisat3n_dna_split_reads_summary_R2_tar,
            dedup_stats = Merge_sort_analyze.dedup_stats_tar,
            chromatin_contact_stats = Merge_sort_analyze.chromatin_contact_stats,
            allc_uniq_reads_stats = Merge_sort_analyze.allc_uniq_reads_stats,
            unique_reads_cgn_extraction_tbi = Merge_sort_analyze.extract_allc_output_tbi_tar,
            plate_id = plate_id,
            docker = docker_prefix + m3c_yap_hisat_docker,
            cromwell_root_dir = cromwell_root_dir,
            cloud_provider = cloud_provider
    }

    meta {
        allowNestedInputs: true
    }

    output {
        File MappingSummary = Summary.mapping_summary
        Array[File] reference_version = Hisat_paired_end.reference_version
        Array[File] name_sorted_bam_array = Summary_PerCellOutput.name_sorted_bam_array
        Array[File] unique_reads_cgn_extraction_allc_array = Summary_PerCellOutput.unique_reads_cgn_extraction_allc_array
        Array[File] unique_reads_cgn_extraction_tbi_array = Summary_PerCellOutput.unique_reads_cgn_extraction_tbi_array
        Array[File] all_reads_3C_contacts_array = Summary_PerCellOutput.all_reads_3C_contacts_array
        Array[File] unique_reads_cgn_extraction_allc_extract_array = Summary_PerCellOutput.unique_reads_cgn_extraction_allc_extract_array
        Array[File] unique_reads_cgn_extraction_tbi_extract_array = Summary_PerCellOutput.unique_reads_cgn_extraction_tbi_extract_array
    }
}

task Demultiplexing {
  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_id
    Int batch_number
    Int min_threshold = 100
    Int max_threshold = 6000000
    String docker

    Int disk_size = 1000
    Int mem_size = 10
    Int preemptible_tries = 2
    Int cpu = 8
  }

  command <<<
    set -euo pipefail
    WORKING_DIR=`pwd`

    # Cat files for each r1, r2
    cat ~{sep=' ' fastq_input_read1} > $WORKING_DIR/r1.fastq.gz
    cat ~{sep=' ' fastq_input_read2} > $WORKING_DIR/r2.fastq.gz

    # Run cutadapt
    /opt/conda/bin/cutadapt -Z -e 0.01 --no-indels -j 8 \
    -g file:~{random_primer_indexes} \
    -o ~{plate_id}-{name}-R1.fq.gz \
    -p ~{plate_id}-{name}-R2.fq.gz \
    $WORKING_DIR/r1.fastq.gz \
    $WORKING_DIR/r2.fastq.gz \
    > $WORKING_DIR/~{plate_id}.stats.txt

    # Remove the fastq files that end in unknown-R1.fq.gz and unknown-R2.fq.gz
    rm $WORKING_DIR/*-unknown-R{1,2}.fq.gz

    python3 <<CODE
    import re
    import os

    # Parsing stats.txt file
    working_dir = os.getcwd()
    stats_file_path = os.path.join(working_dir, '~{plate_id}.stats.txt')
    adapter_counts = {}
    with open(stats_file_path, 'r') as file:
        content = file.read()

    adapter_matches = re.findall(r'=== First read: Adapter (\w+) ===\n\nSequence: .+; Type: .+; Length: \d+; Trimmed: (\d+) times', content)
    for adapter_match in adapter_matches:
        adapter_name = adapter_match[0]
        trimmed_count = int(adapter_match[1])
        adapter_counts[adapter_name] = trimmed_count

    # Removing fastq files with trimmed reads greater than 10000000 or less than 100
    for filename in os.listdir(working_dir):
        if filename.endswith('.fq.gz'):
            file_path = os.path.join(working_dir, filename)
            adapter_name = re.search(r'([A-Za-z]\d+)-R', filename).group(1)
            if adapter_name:
                if adapter_name in adapter_counts:
                    if adapter_counts[adapter_name] < ~{min_threshold} or adapter_counts[adapter_name] > ~{max_threshold}:
                        print("Removing ", file_path, " with count equal to ", adapter_counts[adapter_name])
                        os.remove(file_path)
    CODE
    
    # Check if the number of *R1.fq.gz files is 0
    if [[ $(ls | grep "\-R1.fq.gz" | wc -l) -eq 0 ]]; then
        echo "Error: No files found. All fastq files were removed. Exiting."
        exit 1
    fi

    # Batch the fastq files into folders of batch_number size
    R1_files=($(ls $WORKING_DIR | grep "\-R1.fq.gz"))
    R2_files=($(ls $WORKING_DIR | grep "\-R2.fq.gz"))
    batch_number=~{batch_number}
    total_files=${#R1_files[@]}
    echo "Total files: $total_files"

    if [[ $total_files -lt $batch_number ]]; then
        echo "Warning: Number of files is less than the batch number. Updating batch number to $total_files."
        batch_number=$total_files
    fi

    for i in $(seq 1 "${batch_number}"); do  # Use seq for reliable brace expansion
        mkdir -p "batch${i}"  # Combine batch and i, use -p to create parent dirs
    done

    # Counter for the folder index and create emptycells file
    folder_index=1

    # Distribute the FASTQ files and create TAR files
    for file in "${R1_files[@]}"; do
        sample_id=$(basename "$file" "-R1.fq.gz")
        r2_file="${sample_id}-R2.fq.gz"
         
        mv $WORKING_DIR/$file batch$((folder_index))/$file
        mv $WORKING_DIR/$r2_file batch$((folder_index))/$r2_file
        # Increment the counter
        folder_index=$(( (folder_index % $batch_number) + 1 ))
    done

    # Tar up files per batch
    echo "TAR files"
    for i in $(seq 1 "${batch_number}"); do
        tar -cf - $WORKING_DIR/batch${i}/*.fq.gz | pigz > ~{plate_id}.${i}.cutadapt_output_files.tar.gz
    done
  >>>
  
  runtime {
    docker: docker
    disks: "local-disk ${disk_size} SSD"
    cpu: cpu
    memory: "${mem_size} GiB"
    preemptible: preemptible_tries
  }

  output {
    Array[File] tarred_demultiplexed_fastqs = glob("*.tar.gz")
    File stats = "~{plate_id}.stats.txt"
    }
}

task Hisat_paired_end {
    input {
        File tarred_demultiplexed_fastqs
        File tarred_index_files
        File genome_fa
        File chromosome_sizes
        String plate_id
        String docker
        String cromwell_root_dir
        String cloud_provider

        String r1_adapter
        String r2_adapter
        Int r1_left_cut
        Int r1_right_cut
        Int r2_left_cut
        Int r2_right_cut
        Int min_read_length
        Int disk_size = 1000
        Int cpu = 48
        Int mem_size = 64
        Int preemptible_tries = 2
        String cpu_platform =  "Intel Ice Lake"
    }

    command <<<
        set -euo pipefail
        set -x
        lscpu
        WORKING_DIR=`pwd`

        # check genomic reference version and print to output txt file
        STRING=~{genome_fa}
        BASE=$(basename $STRING .fa)

        echo "The reference is $BASE" > ~{plate_id}.reference_version.txt

        # untar the index files for hisat task
        start=$(date +%s)
        echo "Untarring tarred_index_files"
        pigz -dc ~{tarred_index_files} | tar -xf -
        rm ~{tarred_index_files}
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to untar tarred_index_files: $elapsed seconds"

        # get the basename of the genome_fa file
        cp ~{genome_fa} .
        genome_fa_basename=$(basename ~{genome_fa} .fa)

        start=$(date +%s)
        echo "samtools faidx $genome_fa_basename.fa"
        samtools faidx $genome_fa_basename.fa
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to samtools faidx: $elapsed seconds"

        min_read_length=~{min_read_length}

        # untar the demultiplexed fastqs for sort and trim task
        start=$(date +%s)
        echo "Untar demultiplexed fastqs"
        pigz -dc ~{tarred_demultiplexed_fastqs} | tar -xf -
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to untar: $elapsed seconds"

        if [ ~{cloud_provider} = "gcp" ]; then
            batch_dir="~{cromwell_root_dir}~{cromwell_root_dir}/batch*/"
        else
            batch_dir="~{cromwell_root_dir}/*/*/*/*/*~{cromwell_root_dir}/*/*/*/*/batch*/"
        fi

        task() {
          local file=$1
          sample_id=$(basename "$file" "-R1.fq.gz")
          echo $sample_id

          r2_file="${sample_id}-R2.fq.gz"
          r1_file="${sample_id}-R1.fq.gz"
          cp $batch_dir/"$r1_file" .
          cp $batch_dir/"$r2_file" .

          # sort
          start=$(date +%s)
          echo "Run sort r1"
          zcat "$r1_file" | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > "${sample_id}-R1_sorted.fq"
          end=$(date +%s)
          elapsed=$((end - start))
          echo "Elapsed time to run sort r1: $elapsed seconds"

          # sort
          start=$(date +%s)
          echo "Run sort r2"
          zcat "$r2_file" | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > "${sample_id}-R2_sorted.fq"
          end=$(date +%s)
          elapsed=$((end - start))
          echo "Elapsed time to run sort r2: $elapsed seconds"

          # trim using cutadapt
          start=$(date +%s)
          echo "Run cutadapt"
          /opt/conda/bin/cutadapt \
          -a R1Adapter=~{r1_adapter} \
          -A R2Adapter=~{r2_adapter} \
          --report=minimal -O 6 -q 20 \
          -u ~{r1_left_cut} \
          -u -~{r1_right_cut} \
          -U ~{r2_left_cut} \
          -U -~{r2_right_cut} \
          -Z \
          -m ${min_read_length}:${min_read_length} \
          --pair-filter 'both' \
          -o ${sample_id}-R1_trimmed.fq.gz \
          -p ${sample_id}-R2_trimmed.fq.gz \
          ${sample_id}-R1_sorted.fq ${sample_id}-R2_sorted.fq \
          > ${sample_id}.trimmed.stats.txt
          end=$(date +%s)
          elapsed=$((end - start))
          echo "Elapsed time to run cutadapt: $elapsed seconds"

          # hisat run
          start=$(date +%s)
          echo "Run hisat"
          if [ ~{cloud_provider} = "gcp" ]; then
            hisat_index_file_dir="~{cromwell_root_dir}/$genome_fa_basename"
          else
            hisat_index_file_dir="$WORKING_DIR/$genome_fa_basename"
          fi

          hisat-3n $hisat_index_file_dir \
          -q \
          -1 ${sample_id}-R1_trimmed.fq.gz \
          -2 ${sample_id}-R2_trimmed.fq.gz \
          --directional-mapping-reverse --base-change C,T \
          --no-repeat-index \
          --no-spliced-alignment \
          --no-temp-splicesite \
          -t \
          --new-summary \
          --summary-file ${sample_id}.hisat3n_dna_summary.txt \
          --threads 8 | samtools view -b -q 0 -o "${sample_id}.hisat3n_dna.unsort.bam"
          end=$(date +%s)
          elapsed=$((end - start))
          echo "Elapsed time to run hisat: $elapsed seconds"

          # call separate_unique_and_multi_align_reads
          start=$(date +%s)
          echo "Run separate_unique_and_multi_align_reads"
          python3 -c 'from cemba_data.hisat3n import separate_unique_and_multi_align_reads;separate_unique_and_multi_align_reads(in_bam_path="'"$sample_id"'.hisat3n_dna.unsort.bam", out_unique_path="'"$sample_id"'.hisat3n_dna.unique_aligned.bam", out_multi_path="'"$sample_id"'.hisat3n_dna.multi_aligned.bam", out_unmappable_path="'"$sample_id"'.hisat3n_dna.unmapped.fastq", unmappable_format="fastq", mapq_cutoff=10, qlen_cutoff='"$min_read_length"')'
          end=$(date +%s)
          elapsed=$((end - start))
          echo "Elapsed time to run separate_unique_and_multi_align_reads: $elapsed seconds"

          # call split_hisat3n_unmapped_reads
          start=$(date +%s)
          echo "Run split_hisat3n_unmapped_reads"
          python3 -c 'from cemba_data.hisat3n import *;split_hisat3n_unmapped_reads(fastq_path="'"$sample_id"'.hisat3n_dna.unmapped.fastq",output_prefix="'"$sample_id"'.hisat3n_dna.split_reads",min_length='"$min_read_length"')'
          end=$(date +%s)
          elapsed=$((end - start))
          echo "Elapsed time to run split_hisat3n_unmapped_reads: $elapsed seconds"

          rm ${sample_id}-R1.fq.gz ${sample_id}-R2.fq.gz
          rm ${sample_id}-R1_sorted.fq ${sample_id}-R2_sorted.fq
          rm ${sample_id}-R1_trimmed.fq.gz ${sample_id}-R2_trimmed.fq.gz
          rm ${sample_id}.hisat3n_dna.unsort.bam ${sample_id}.hisat3n_dna.multi_aligned.bam
          rm ${sample_id}.hisat3n_dna.unmapped.fastq
       }


      R1_files=($(ls $batch_dir | grep "\-R1.fq.gz"))
      R2_files=($(ls $batch_dir | grep "\-R2.fq.gz"))

      # run 6 instances of task in parallel
      for file in "${R1_files[@]}"; do
        (
          echo "starting task $file.."
          task "$file"
          sleep $(( (RANDOM % 3) + 1))
        ) &
        # allow to execute up to 6 jobs in parallel
        if [[ $(jobs -r -p | wc -l) -ge 6 ]]; then
          wait -n
        fi
      done

      wait
      echo "Tasks all done."
      du -h *

      ####################################
      ## make sure that the number of output bams equals the length of R1_files
      # Count the number of *.hisat3n_dna.unique_aligned.bam files
      bam_count=$(find . -maxdepth 1 -type f -name '*.hisat3n_dna.unique_aligned.bam' | wc -l)
      fastq_counts=$(find . -maxdepth 1 -type f -name '*.split_reads*.fastq' | wc -l)

      # Get the length of the array ${R1_files[@]}
      array_length=${#R1_files[@]}

      # Check if the count of *.hisat3n_dna.unique_aligned.bam files matches the length of the array ${R1_files[@]}
      if [ "$bam_count" -ne "$array_length" ]; then
         echo "Error: Number of BAM files does not match the length of the array."
         exit 1
      fi

      # Check if the count of FASTQ files matches the length of the array ${R1_files[@]}
      if [ "$fastq_counts" -ne  "$((2 * array_length))" ]; then
         echo "Error: Number of FASTQ files: $fastq_count does not match the 2 * length of the array: ${#R1_files[@]}."
         exit 1
      fi

      echo "Number of BAM and FASTQ files matches the length of the array."
      ####################################

      # tar up stats
      echo "Tar up stats"
      start=$(date +%s)
      tar -cf - *.trimmed.stats.txt | pigz > ~{plate_id}.trimmed_stats_files.tar.gz
      tar -cf - *.hisat3n_dna_summary.txt | pigz > ~{plate_id}.hisat3n_paired_end_stats_files.tar.gz
      end=$(date +%s)
      elapsed=$((end - start))
      echo "Elapsed time to run tar stats $elapsed seconds"

      # tar up the uniqe bams
      echo "Tar up unique bams"
      start=$(date +%s)
      tar -cf - *.hisat3n_dna.unique_aligned.bam | pigz > ~{plate_id}.hisat3n_paired_end_unique_bam_files.tar.gz
      end=$(date +%s)
      elapsed=$((end - start))
      echo "Elapsed time to run tar unique bams $elapsed seconds"

      # tar up the split fastq files
      echo "Tar up fastqs"
      start=$(date +%s)
      tar -cf - *.split_reads*.fastq | pigz > ~{plate_id}.hisat3n_paired_end_split_fastq_files.tar.gz
      end=$(date +%s)
      elapsed=$((end - start))
      echo "Elapsed time to run tar fastqs $elapsed seconds"

    >>>

    runtime {
        docker: docker
        disks: "local-disk ${disk_size} SSD"
        cpu: cpu
        memory: "${mem_size} GiB"
        cpuPlatform: cpu_platform
        preemptible: preemptible_tries
    }

    output {
        File trim_stats_tar = "~{plate_id}.trimmed_stats_files.tar.gz"
        File hisat3n_paired_end_stats_tar = "~{plate_id}.hisat3n_paired_end_stats_files.tar.gz"
        File unique_bam_tar = "~{plate_id}.hisat3n_paired_end_unique_bam_files.tar.gz"
        File split_fq_tar = "~{plate_id}.hisat3n_paired_end_split_fastq_files.tar.gz"
        File reference_version = "~{plate_id}.reference_version.txt"
    }
}

task Hisat_single_end {
    input {
        File split_fq_tar
        File genome_fa
        File tarred_index_files
        String plate_id
        String docker
        String cromwell_root_dir
        String cloud_provider

        Int disk_size = 1000
        Int mem_size = 64
        Int cpu = 32
        Int preemptible_tries = 2
        String cpu_platform =  "Intel Ice Lake"
    }

    command <<<
        set -euo pipefail
        set -x
        lscpu
        WORKING_DIR=`pwd`

        # untar the tarred index files
        echo "Untar tarred_index_files"
        start=$(date +%s)
        pigz -dc ~{tarred_index_files} | tar -xf -
        rm ~{tarred_index_files}
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to untar tarred_index_files: $elapsed seconds"

        cp ~{genome_fa} .

        #get the basename of the genome_fa file
        echo "samtools faidx"
        start=$(date +%s)
        genome_fa_basename=$(basename ~{genome_fa} .fa)
        samtools faidx $genome_fa_basename.fa
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to samtools faidx: $elapsed seconds"

        # untar the unmapped fastq files
        echo "Untar split_fq_tar"
        start=$(date +%s)
        pigz -dc ~{split_fq_tar} | tar -xf -
        rm ~{split_fq_tar}
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to untar split_fq_tar: $elapsed seconds"

        # make directories
        mkdir -p ~{cromwell_root_dir}/merged_sort_bams
        mkdir -p ~{cromwell_root_dir}/read_overlap

        # define lists of r1 and r2 fq files
        R1_files=($(ls | grep "\.hisat3n_dna.split_reads.R1.fastq"))
        R2_files=($(ls | grep "\.hisat3n_dna.split_reads.R2.fastq"))

        task() {
          BASE=$(basename "$file" ".hisat3n_dna.split_reads.R1.fastq")
          echo $BASE
          echo "Running hisat on sample_id_R1" $BASE

          echo "Hisat 3n R1"
          start=$(date +%s)

          if [ ~{cloud_provider} = "gcp" ]; then
            hisat_index_file_dir="~{cromwell_root_dir}/$genome_fa_basename"
          else
            hisat_index_file_dir="$WORKING_DIR/$genome_fa_basename"
          fi


          # hisat on R1 single end
          hisat-3n $hisat_index_file_dir \
          -q \
          -U ${BASE}.hisat3n_dna.split_reads.R1.fastq \
          -S ${BASE}.hisat3n_dna.split_reads.R1.sam --directional-mapping-reverse --base-change C,T \
          --no-repeat-index \
          --no-spliced-alignment \
          --no-temp-splicesite \
          -t \
          --new-summary \
          --summary-file ${BASE}.hisat3n_dna_split_reads_summary.R1.txt \
          --threads 8

         end=$(date +%s)
         elapsed=$((end - start))
         echo "Elapsed time to run $elapsed seconds"
         echo "Finish running hisat on sample_id_R1" $BASE

         echo "Hisat 3n R2"
         start=$(date +%s)
         echo "Running hisat on sample_id_R2" $BASE

         # hisat on R2 single end
         hisat-3n $hisat_index_file_dir \
         -q \
         -U ${BASE}.hisat3n_dna.split_reads.R2.fastq \
         -S ${BASE}.hisat3n_dna.split_reads.R2.sam --directional-mapping --base-change C,T \
         --no-repeat-index \
         --no-spliced-alignment \
         --no-temp-splicesite \
         -t --new-summary \
         --summary-file ${BASE}.hisat3n_dna_split_reads_summary.R2.txt \
         --threads 8

         end=$(date +%s)
         elapsed=$((end - start))
         echo "Elapsed time to run $elapsed seconds"
         echo "Finish running hisat on sample_id_R2" $BASE

         # samtools merge
         echo "samtools merge R1 and R2"
         start=$(date +%s)
         samtools merge -o ${BASE}.name_merged.sam ${BASE}.hisat3n_dna.split_reads.R1.sam ${BASE}.hisat3n_dna.split_reads.R2.sam -@8
         end=$(date +%s)
         elapsed=$((end - start))
         echo "Elapsed time to run samtools merge $elapsed seconds"

         # samtools sort
         echo "samtools sort R1 and R2"
         start=$(date +%s)
         samtools sort -n -@8 -m1g ${BASE}.name_merged.sam -o ${BASE}.name_sorted.bam
         end=$(date +%s)
         elapsed=$((end - start))
         echo "Elapsed time to run samtools sort $elapsed seconds"

         # samtools filter bam
         echo "samtools -q 10"
         start=$(date +%s)
         samtools view -q 10 ${BASE}.name_sorted.bam -o ${BASE}.name_sorted.filtered.bam
         end=$(date +%s)
         elapsed=$((end - start))
         echo "Elapsed time to run samtools -q 10 $elapsed seconds"

         if [ ~{cloud_provider} = "gcp" ]; then
            bam_path_prefix="~{cromwell_root_dir}"
         else
            bam_path_prefix=$WORKING_DIR
         fi

         # remove_overlap_read_parts
         echo "call remove_overlap_read_parts"
         start=$(date +%s)
         python3 -c 'from cemba_data.hisat3n import *;import os;remove_overlap_read_parts(in_bam_path="'"$BASE"'.name_sorted.filtered.bam",out_bam_path="'"$BASE"'.hisat3n_dna.split_reads.read_overlap.bam")'
         end=$(date +%s)
         elapsed=$((end - start))
         echo "Elapsed time to run remove overlap $elapsed seconds"

         # remove files
         rm ${BASE}.hisat3n_dna.split_reads.R1.fastq ${BASE}.hisat3n_dna.split_reads.R2.fastq
         rm ${BASE}.hisat3n_dna.split_reads.R1.sam ${BASE}.hisat3n_dna.split_reads.R2.sam
         rm ${BASE}.name_merged.sam
         rm ${BASE}.name_sorted.bam
         rm ${BASE}.name_sorted.filtered.bam
        }

        # run 4 instances in parallel each with 8 threads
        for file in "${R1_files[@]}"; do
          (
            echo "starting task $file.."
            task "$file"
            sleep $(( (RANDOM % 3) + 1))
          ) &
          if [[ $(jobs -r -p | wc -l) -ge 4 ]]; then
            wait -n
          fi
        done

        wait
        echo "All done running tasks."
        du -h *

        ####################################
        ## make sure that the number of output bams equals the length of R1_files
        # Count the number of bam files
        bam_count=$(find . -maxdepth 1 -type f -name '*read_overlap.bam' | wc -l)

        # Get the length of the array ${R1_files[@]}
        array_length=${#R1_files[@]}

        # Check if the count of bams matches the length of the array ${R1_files[@]}
        if [ "$bam_count" -ne "$array_length" ]; then
           echo "Error: Number of BAM files does not match the length of the array."
           exit 1
        fi
        echo "Number of BAM files matches the length of the array."
        ####################################

        echo "Tar up summary text files"
        start=$(date +%s)
        # tar up the r1 and r2 stats files -p to set number of threads
        tar -cf - *.hisat3n_dna_split_reads_summary.R1.txt | pigz > ~{plate_id}.hisat3n_dna_split_reads_summary.R1.tar.gz
        tar -cf - *.hisat3n_dna_split_reads_summary.R2.txt | pigz > ~{plate_id}.hisat3n_dna_split_reads_summary.R2.tar.gz
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to run tar summary text files $elapsed seconds"

        # tar up read overlap files
        echo "Tar up read_overlap bams"
        start=$(date +%s)
        tar -cf - *read_overlap.bam | pigz > ~{plate_id}.remove_overlap_read_parts.tar.gz
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to tar read_overlap bams $elapsed seconds"
    >>>

    runtime {
        docker: docker
        disks: "local-disk ${disk_size} SSD"
        cpu: cpu
        memory: "${mem_size} GiB"
        cpuPlatform: cpu_platform
        preemptible: preemptible_tries
    }

    output {
         File hisat3n_dna_split_reads_summary_R1_tar = "~{plate_id}.hisat3n_dna_split_reads_summary.R1.tar.gz"
         File hisat3n_dna_split_reads_summary_R2_tar = "~{plate_id}.hisat3n_dna_split_reads_summary.R2.tar.gz"
         File remove_overlaps_output_bam_tar = "~{plate_id}.remove_overlap_read_parts.tar.gz"

    }
}

task Merge_sort_analyze {
    input {
        String plate_id
        File paired_end_unique_tar
        File read_overlap_tar
        String docker
        String cromwell_root_dir
        String cloud_provider

        #input for allcools bam-to-allc
        File genome_fa
        String genome_base = basename(genome_fa)
        Int num_upstr_bases
        Int num_downstr_bases
        Int compress_level
        File chromosome_sizes

        String cpu_platform = "Intel Ice Lake"
        Int disk_size = 1000
        Int mem_size = 40
        Int cpu = 16
        Int preemptible_tries = 2
    }

    command <<<
      set -euo pipefail
      set -x
      lscpu

      WORKING_DIR=`pwd`

      # unzip tars
      echo "Untar paired_end_unique_tar"
      start=$(date +%s)
      pigz -dc ~{paired_end_unique_tar} | tar -xf -
      rm ~{paired_end_unique_tar}
      end=$(date +%s)
      elapsed=$((end - start))
      echo "Elapsed time to untar paired_end_unique_tar: $elapsed seconds"

      echo "Untar read_overlap_tar"
      start=$(date +%s)
      pigz -dc ~{read_overlap_tar} | tar -xf -
      rm ~{read_overlap_tar}
      end=$(date +%s)
      elapsed=$((end - start))
      echo "Elapsed time to untar read_overlap_tar: $elapsed seconds"

      # reference and index
      start=$(date +%s)
      echo "Reference and index fasta"
      mkdir reference
      cp ~{genome_fa} reference
      ls reference
      samtools faidx reference/*.fa
      end=$(date +%s)
      elapsed=$((end - start))
      echo "Elapsed time to index fasta $elapsed seconds"

      # define lists of r1 and r2 fq files
      UNIQUE_BAMS=($(ls | grep "\.hisat3n_dna.unique_aligned.bam"))
      SPLIT_BAMS=($(ls | grep "\.hisat3n_dna.split_reads.read_overlap.bam"))

      # for allcools bam-to-allc
      if [ ~{num_upstr_bases} -eq 0 ]; then
         mcg_context=CGN
      else
         mcg_context=HCGN
      fi

      # make directories
      mkdir ~{cromwell_root_dir}/output_bams
      mkdir ~{cromwell_root_dir}/temp
      mkdir ~{cromwell_root_dir}/allc-${mcg_context}

      task() {
        local file=$1
        sample_id=$(basename "$file" ".hisat3n_dna.unique_aligned.bam")
        echo $sample_id

        start=$(date +%s)
        echo "Merge all unique_aligned and read_overlap"
        samtools merge -f "${sample_id}.hisat3n_dna.all_reads.bam" "${sample_id}.hisat3n_dna.unique_aligned.bam" "${sample_id}.hisat3n_dna.split_reads.read_overlap.bam" -@4
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to run merge $elapsed seconds"

        start=$(date +%s)
        echo "Sort all reads by name"
        samtools sort -n -@4 -m1g -o "${sample_id}.hisat3n_dna.all_reads.name_sort.bam" "${sample_id}.hisat3n_dna.all_reads.bam"
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to run sort by name $elapsed seconds"

        start=$(date +%s)
        echo "Sort all reads by position"
        samtools sort -O BAM -@4 -m1g -o "${sample_id}.hisat3n_dna.all_reads.pos_sort.bam" "${sample_id}.hisat3n_dna.all_reads.name_sort.bam"
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to run sort by pos $elapsed seconds"

        start=$(date +%s)
        echo "Call Picard remove duplicates"
        name=${sample_id}.hisat3n_dna.all_reads.deduped
        picard MarkDuplicates I=${sample_id}.hisat3n_dna.all_reads.pos_sort.bam O=~{cromwell_root_dir}/output_bams/${name}.bam \
        M=~{cromwell_root_dir}/output_bams/${name}.matrix.txt \
        REMOVE_DUPLICATES=true TMP_DIR=~{cromwell_root_dir}/temp
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to run picard $elapsed seconds"

        start=$(date +%s)
        echo "Call samtools index"
        samtools index ~{cromwell_root_dir}/output_bams/${name}.bam
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to samtools index $elapsed seconds"

        start=$(date +%s)
        echo "Call chromatin contacts from name sorted bams"
        python3 -c 'from cemba_data.hisat3n import *;import os;import glob;call_chromatin_contacts(bam_path="'"$sample_id"'.hisat3n_dna.all_reads.name_sort.bam",contact_prefix="'"$sample_id"'.hisat3n_dna.all_reads",save_raw=False,save_hic_format=True)'
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to chromatin contacts $elapsed seconds"

        if [ ~{cloud_provider} = "gcp" ]; then
            reference_fasta="~{cromwell_root_dir}/reference/~{genome_base}"
          else
            reference_fasta="$WORKING_DIR/reference/~{genome_base}"
        fi

        start=$(date +%s)
        echo "Call allcools bam-to-allc from deduped.bams"
        /opt/conda/bin/allcools bam-to-allc \
        --bam_path ~{cromwell_root_dir}/output_bams/${name}.bam \
        --reference_fasta $reference_fasta \
        --output_path "${sample_id}.allc.tsv.gz" \
        --num_upstr_bases ~{num_upstr_bases} \
        --num_downstr_bases ~{num_downstr_bases} \
        --compress_level ~{compress_level} \
        --save_count_df \
        --convert_bam_strandness
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to allcools bam-to-allc $elapsed seconds"

        start=$(date +%s)
        echo "Call allcools extract-all"
        allcools extract-allc --strandness merge \
        --allc_path ${sample_id}.allc.tsv.gz \
        --output_prefix ~{cromwell_root_dir}/allc-${mcg_context}/${sample_id} \
        --mc_contexts ${mcg_context} \
        --chrom_size_path ~{chromosome_sizes}
        end=$(date +%s)
        elapsed=$((end - start))
        echo "Elapsed time to allcools extract-all $elapsed seconds"

        echo "Remove some bams"
        rm ${sample_id}.hisat3n_dna.all_reads.bam
        rm ${sample_id}.hisat3n_dna.all_reads.pos_sort.bam
        rm ~{cromwell_root_dir}/${sample_id}.hisat3n_dna.split_reads.read_overlap.bam
        rm ~{cromwell_root_dir}/${sample_id}.hisat3n_dna.unique_aligned.bam
      }

      # run 4 instances of task in parallel
      for file in "${UNIQUE_BAMS[@]}"; do
        (
          echo "starting task $file.."
          task "$file"
          sleep $(( (RANDOM % 3) + 1))
        ) &
        # allow to execute up to 4 jobs in parallel
        if [[ $(jobs -r -p | wc -l) -ge 4 ]]; then
          wait -n
        fi
      done

      wait
      echo "Tasks all done."
      du -h *

      ####################################
      ## make sure that the number of output bams equals the length of UNIQUE_BAMS
      # Count the number of *.hisat3n_dna.unique_aligned.bam files
      bam_count=$(find . -maxdepth 1 -type f -name '*.hisat3n_dna.all_reads.name_sort.bam' | wc -l)
      contact_count=$(find . -maxdepth 1 -type f -name '*.hisat3n_dna.all_reads.3C.contact.tsv.gz' | wc -l)

      # Get the length of the array ${UNIQUE_BAMS[@]}
      array_length=${#UNIQUE_BAMS[@]}

      # Check if the count of bams matches the length of the array ${UNIQUE_BAMS[@]}
      if [ "$bam_count" -ne "$array_length" ]; then
          echo "Error: Number of BAM files does not match the length of the array."
          exit 1
      fi
      # Check if the count of tsv files matches the length of the array ${UNIQUE_BAMS[@]}
      if [ "$contact_count" -ne "$array_length" ]; then
          echo "Error: Number of hisat3n_dna.all_reads.3C.contact.tsv.gz files does not match the length of the array."
          exit 1
      fi
      echo "Number of output files matches the length of the array."
      ####################################

      echo "recursively ls'sing cromwell root again"
      ls -lR ~{cromwell_root_dir}

      echo "Tar files."
      tar -cf - ~{cromwell_root_dir}/output_bams/*.matrix.txt | pigz > ~{plate_id}.dedup_unique_bam_and_index_unique_bam_stats.tar.gz
      tar -cf - *.hisat3n_dna.all_reads.name_sort.bam | pigz > ~{plate_id}.hisat3n_dna.all_reads.name_sort.tar.gz

      # tar outputs of call_chromatin_contacts
      tar -cf - *.hisat3n_dna.all_reads.3C.contact.tsv.gz | pigz > ~{plate_id}.hisat3n_dna.all_reads.3C.contact.tar.gz
      tar -cf - *.hisat3n_dna.all_reads.dedup_contacts.tsv.gz | pigz > ~{plate_id}.hisat3n_dna.all_reads.dedup_contacts.tar.gz
      tar -cf - *.hisat3n_dna.all_reads.contact_stats.csv | pigz > ~{plate_id}.chromatin_contact_stats.tar.gz

      # tar outputs of allcools
      tar -cf - *.allc.tsv.gz | pigz > ~{plate_id}.allc.tsv.tar.gz
      tar -cf - *.allc.tsv.gz.tbi | pigz > ~{plate_id}.allc.tbi.tar.gz
      tar -cf -  *.allc.tsv.gz.count.csv | pigz > ~{plate_id}.allc.count.tar.gz
      tar -cf -  ~{cromwell_root_dir}/allc-${mcg_context}/*.gz | pigz > ~{plate_id}.extract-allc.tar.gz
      tar -cf -  ~{cromwell_root_dir}/allc-${mcg_context}/*.tbi | pigz > ~{plate_id}.extract-allc_tbi.tar.gz
    >>>

    runtime {
        docker: docker
        disks: "local-disk ${disk_size} SSD"
        cpu: cpu
        memory: "${mem_size} GiB"
        cpuPlatform: cpu_platform
        preemptible: preemptible_tries
    }

     output {
        File allc = "~{plate_id}.allc.tsv.tar.gz"
        File tbi = "~{plate_id}.allc.tbi.tar.gz"
        File all_reads_dedup_contacts = "~{plate_id}.hisat3n_dna.all_reads.dedup_contacts.tar.gz"
        File all_reads_3C_contacts = "~{plate_id}.hisat3n_dna.all_reads.3C.contact.tar.gz"
        File name_sorted_bam = "~{plate_id}.hisat3n_dna.all_reads.name_sort.tar.gz"
        File dedup_stats_tar = "~{plate_id}.dedup_unique_bam_and_index_unique_bam_stats.tar.gz"
        File chromatin_contact_stats = "~{plate_id}.chromatin_contact_stats.tar.gz"
        File allc_uniq_reads_stats = "~{plate_id}.allc.count.tar.gz"
        File extract_allc_output_tbi_tar = "~{plate_id}.extract-allc_tbi.tar.gz"
        File extract_allc_output_allc_tar  = "~{plate_id}.extract-allc.tar.gz"
     }
}

task Summary_PerCellOutput {
    input {
        # to generate list of files per file type
        Array[File] name_sorted_bams
        Array[File] unique_reads_cgn_extraction_allc
        Array[File] unique_reads_cgn_extraction_tbi
        Array[File] all_reads_3C_contacts
        Array[File] unique_reads_cgn_extraction_allc_extract
        Array[File] unique_reads_cgn_extraction_tbi_extract
        String cromwell_root_dir

        String docker
        String plate_id
        Int disk_size = 1000
        Int mem_size = 30
        Int cpu = 16
    }

    command <<<
        set -euo pipefail
        set -x

        # Set root_dir to current working directory
        root_dir=$(pwd)
        echo "This is the root directory " $root_dir

        extract_and_remove() {
            if [ $# -eq 0 ];
                then
                    echo "No files exist"
                    return
            fi

            for tarred_file in "${@}"; do
                dir_name=`basename "${tarred_file%.tar.gz}"`
                echo $dir_name
                mkdir -p "$root_dir"/"$dir_name"
                pigz -dc "$tarred_file" | tar -xvf - -C "$root_dir"/"$dir_name"
                rm "$tarred_file"
            done
        }

        # output files at a cell level
        echo "Untar files needed at per cell level"
        extract_and_remove ~{sep=' ' name_sorted_bams}
        extract_and_remove ~{sep=' ' unique_reads_cgn_extraction_allc}
        extract_and_remove ~{sep=' ' unique_reads_cgn_extraction_tbi}
        extract_and_remove ~{sep=' ' all_reads_3C_contacts}
        extract_and_remove ~{sep=' ' unique_reads_cgn_extraction_allc_extract}
        extract_and_remove ~{sep=' ' unique_reads_cgn_extraction_tbi_extract}

    >>>

    runtime {
        docker: docker
        disks: "local-disk ${disk_size} SSD"
        cpu: cpu
        memory: "${mem_size} GiB"
    }

    output {
        Array[File] name_sorted_bam_array = glob("~{plate_id}.hisat3n_dna.all_reads.name_sort/*")
        Array[File] unique_reads_cgn_extraction_allc_array = glob("~{plate_id}.allc.tsv/*")
        Array[File] unique_reads_cgn_extraction_tbi_array = glob("~{plate_id}.allc.tbi/*")
        Array[File] all_reads_3C_contacts_array = glob("~{plate_id}.hisat3n_dna.all_reads.3C.contact/*")
        Array[File] unique_reads_cgn_extraction_allc_extract_array = glob("~{plate_id}.extract-allc/~{cromwell_root_dir}/allc-CGN/*")
        Array[File] unique_reads_cgn_extraction_tbi_extract_array = glob("~{plate_id}.extract-allc_tbi/~{cromwell_root_dir}/allc-CGN/*")
    }
}

task Summary {
    input {
        Array[File] trimmed_stats
        Array[File] hisat3n_stats
        Array[File] r1_hisat3n_stats
        Array[File] r2_hisat3n_stats
        Array[File] dedup_stats
        Array[File] chromatin_contact_stats
        Array[File] allc_uniq_reads_stats
        Array[File] unique_reads_cgn_extraction_tbi
        String plate_id
        String cromwell_root_dir
        String cloud_provider

        String docker
        Int disk_size = 80
        Int mem_size = 5
        Int preemptible_tries = 3
        Int cpu = 4
    }
    command <<<
        set -euo pipefail

        WORKING_DIR=`pwd`

        if [ ~{cloud_provider} = "gcp" ]; then
            base_directory=~{cromwell_root_dir}
            matrix_files_dir="~{cromwell_root_dir}~{cromwell_root_dir}/output_bams"
            allc_index_dir="~{cromwell_root_dir}~{cromwell_root_dir}/allc-*"
        else
            base_directory=$WORKING_DIR
            matrix_files_dir="$WORKING_DIR~{cromwell_root_dir}/output_bams"
            allc_index_dir="$WORKING_DIR~{cromwell_root_dir}/allc-*"
        fi

        mkdir $base_directory/fastq
        mkdir $base_directory/bam
        mkdir $base_directory/allc
        mkdir $base_directory/hic

        extract_and_remove() {
            if [ $# -eq 0 ];
                then
                    echo "No files exist"
                    return
            fi
            for tar in "${@}"; do
                tar -xvf "$tar"
                rm "$tar"
            done
        }

        extract_and_remove ~{sep=' ' trimmed_stats}
        extract_and_remove ~{sep=' ' hisat3n_stats}
        extract_and_remove ~{sep=' ' r1_hisat3n_stats}
        extract_and_remove ~{sep=' ' r2_hisat3n_stats}
        extract_and_remove ~{sep=' ' dedup_stats}
        extract_and_remove ~{sep=' ' chromatin_contact_stats}
        extract_and_remove ~{sep=' ' allc_uniq_reads_stats}
        extract_and_remove ~{sep=' ' unique_reads_cgn_extraction_tbi}

        mv *.trimmed.stats.txt $base_directory/fastq
        mv *.hisat3n_dna_summary.txt *.hisat3n_dna_split_reads_summary.R1.txt *.hisat3n_dna_split_reads_summary.R2.txt $base_directory/bam
        mv $matrix_files_dir/*.hisat3n_dna.all_reads.deduped.matrix.txt $base_directory/bam
        mv *.hisat3n_dna.all_reads.contact_stats.csv $base_directory/hic
        mv *.allc.tsv.gz.count.csv $base_directory/allc
        mv $allc_index_dir/*.allc.tsv.gz.tbi $base_directory/allc

        python3 -c 'from cemba_data.hisat3n import *;snm3c_summary()'
        mv MappingSummary.csv.gz ~{plate_id}_MappingSummary.csv.gz

    >>>
    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
    }
    output {
        File mapping_summary = "~{plate_id}_MappingSummary.csv.gz"
    }
}
