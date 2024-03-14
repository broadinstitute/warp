version 1.0

workflow snM3C {

    input {
        Array[File] fastq_input_read1
        Array[File] fastq_input_read2
        File random_primer_indexes
        String plate_id
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
        Int batch_number
        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:2.3"
    }

    # version of the pipeline
    String pipeline_version = "3.0.0"

    call Demultiplexing {
        input:
            fastq_input_read1 = fastq_input_read1,
            fastq_input_read2 = fastq_input_read2,
            random_primer_indexes = random_primer_indexes,
            plate_id = plate_id,
            docker = docker,
            batch_number = batch_number
    }

    scatter(tar in Demultiplexing.tarred_demultiplexed_fastqs) {
        call hisat_paired_end {
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
                docker = docker
        }

        call hisat_single_end {
            input:
                split_fq_tar = hisat_paired_end.split_fq_tar,
                tarred_index_files = tarred_index_files,
                genome_fa = genome_fa,
                plate_id = plate_id,
                docker = docker
        }

        call merge_sort_analyze {
            input:
               paired_end_unique_tar = hisat_paired_end.unique_bam_tar,
               read_overlap_tar = hisat_single_end.remove_overlaps_output_bam_tar,     
               genome_fa = genome_fa, 
               num_upstr_bases = num_upstr_bases,
               num_downstr_bases = num_downstr_bases,
               compress_level = compress_level,
               chromosome_sizes = chromosome_sizes,
               plate_id = plate_id,
               docker = docker
        }
    }

    call summary {
        input:
            trimmed_stats = hisat_paired_end.trim_stats_tar,
            hisat3n_stats = hisat_paired_end.hisat3n_paired_end_stats_tar,
            r1_hisat3n_stats = hisat_single_end.hisat3n_dna_split_reads_summary_R1_tar,
            r2_hisat3n_stats = hisat_single_end.hisat3n_dna_split_reads_summary_R2_tar,
            dedup_stats = merge_sort_analyze.dedup_stats_tar,
            chromatin_contact_stats = merge_sort_analyze.chromatin_contact_stats,
            allc_uniq_reads_stats = merge_sort_analyze.allc_uniq_reads_stats,
            unique_reads_cgn_extraction_tbi = merge_sort_analyze.extract_allc_output_tbi_tar,
            plate_id = plate_id,
            docker = docker
    }

    output {
        File MappingSummary = summary.mapping_summary
        Array[File] name_sorted_bams = merge_sort_analyze.name_sorted_bam
        Array[File] unique_reads_cgn_extraction_allc= merge_sort_analyze.allc
        Array[File] unique_reads_cgn_extraction_tbi = merge_sort_analyze.tbi
        Array[File] reference_version = hisat_paired_end.reference_version
        Array[File] all_reads_dedup_contacts = merge_sort_analyze.all_reads_dedup_contacts
        Array[File] all_reads_3C_contacts = merge_sort_analyze.all_reads_3C_contacts
        Array[File] chromatin_contact_stats = merge_sort_analyze.chromatin_contact_stats
        Array[File] unique_reads_cgn_extraction_allc_extract = merge_sort_analyze.extract_allc_output_allc_tar
        Array[File] unique_reads_cgn_extraction_tbi_extract = merge_sort_analyze.extract_allc_output_tbi_tar

    }
}

task Demultiplexing {
  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_id
    Int batch_number
    String docker

    Int disk_size = 1000
    Int mem_size = 10
    Int preemptible_tries = 2
    Int cpu = 8
  }

  command <<<
    set -euo pipefail

    # Cat files for each r1, r2
    cat ~{sep=' ' fastq_input_read1} > r1.fastq.gz
    cat ~{sep=' ' fastq_input_read2} > r2.fastq.gz

    /opt/conda/bin/cutadapt -Z -e 0.01 --no-indels -j 8 \
    -g file:~{random_primer_indexes} \
    -o ~{plate_id}-{name}-R1.fq.gz \
    -p ~{plate_id}-{name}-R2.fq.gz \
    r1.fastq.gz \
    r2.fastq.gz \
    > ~{plate_id}.stats.txt

    # remove the fastq files that end in unknown-R1.fq.gz and unknown-R2.fq.gz
    rm *-unknown-R{1,2}.fq.gz

    python3 <<CODE
    import re
    import os

    # Parsing stats.txt file
    stats_file_path = '/cromwell_root/~{plate_id}.stats.txt'
    adapter_counts = {}
    with open(stats_file_path, 'r') as file:
        content = file.read()

    adapter_matches = re.findall(r'=== First read: Adapter (\w+) ===\n\nSequence: .+; Type: .+; Length: \d+; Trimmed: (\d+) times', content)
    for adapter_match in adapter_matches:
        adapter_name = adapter_match[0]
        trimmed_count = int(adapter_match[1])
        adapter_counts[adapter_name] = trimmed_count

    # Removing fastq files with trimmed reads greater than 30
    directory_path = '/cromwell_root'
    threshold = 10000000

    for filename in os.listdir(directory_path):
        if filename.endswith('.fq.gz'):
            file_path = os.path.join(directory_path, filename)
            adapter_name = re.search(r'A(\d+)-R', filename)
            if adapter_name:
                adapter_name = 'A' + adapter_name.group(1)
                if adapter_name in adapter_counts and adapter_counts[adapter_name] > threshold:
                    os.remove(file_path)
                    print(f'Removed file: {filename}')
    CODE

    # Batch the fastq files into folders of batch_number size
    batch_number=~{batch_number}
    for i in $(seq 1 "${batch_number}"); do  # Use seq for reliable brace expansion
        mkdir -p "batch${i}"  # Combine batch and i, use -p to create parent dirs
    done

    # Counter for the folder index
    folder_index=1

    # Define lists of r1 and r2 fq files
    R1_files=($(ls | grep "\-R1.fq.gz"))
    R2_files=($(ls | grep "\-R2.fq.gz"))

    # Distribute the FASTQ files and create TAR files
    for file in "${R1_files[@]}"; do
        sample_id=$(basename "$file" "-R1.fq.gz")
        r2_file="${sample_id}-R2.fq.gz"
        mv $file batch$((folder_index))/$file
        mv $r2_file batch$((folder_index))/$r2_file
        # Increment the counter
        folder_index=$(( (folder_index % $batch_number) + 1 ))
    done
    echo "TAR files"
    for i in $(seq 1 "${batch_number}"); do
        tar -zcvf "~{plate_id}.${i}.cutadapt_output_files.tar.gz" batch${i}/*.fq.gz
    done


    echo "TAR files created successfully."
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

task hisat_paired_end{
    input {
        File tarred_demultiplexed_fastqs
        File tarred_index_files
        File genome_fa
        File chromosome_sizes
        String plate_id
        String docker

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
  
        # check genomic reference version and print to output txt file
        STRING=~{genome_fa}
        BASE=$(basename $STRING .fa)

        echo "The reference is $BASE" > ~{plate_id}.reference_version.txt

        # untar the index files for hisat task
        start=$(date +%s)
        echo "Untarring tarred_index_files"
        tar -zxvf ~{tarred_index_files}
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
        tar -xf ~{tarred_demultiplexed_fastqs}
        end=$(date +%s) 
        elapsed=$((end - start)) 
        echo "Elapsed time to untar: $elapsed seconds"
    
        task() {
          local file=$1
          sample_id=$(basename "$file" "-R1.fq.gz")
          echo $sample_id

          r2_file="${sample_id}-R2.fq.gz"
          r1_file="${sample_id}-R1.fq.gz"
          
          # sort 
          start=$(date +%s)
          echo "Run sort r1"
          zcat /cromwell_root/batch*/"$r1_file" | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > "${sample_id}-R1_sorted.fq"
          end=$(date +%s) 
          elapsed=$((end - start)) 
          echo "Elapsed time to run sort r1: $elapsed seconds"
    
          # sort 
          start=$(date +%s)
          echo "Run sort r1"
          zcat /cromwell_root/batch*/"$r2_file" | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > "${sample_id}-R2_sorted.fq"
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
          hisat-3n /cromwell_root/$genome_fa_basename \
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
          
          rm /cromwell_root/batch*/${sample_id}-R1.fq.gz /cromwell_root/batch*/${sample_id}-R2.fq.gz
          rm ${sample_id}-R1_sorted.fq ${sample_id}-R2_sorted.fq
          rm ${sample_id}-R1_trimmed.fq.gz ${sample_id}-R2_trimmed.fq.gz
          rm ${sample_id}.hisat3n_dna.unsort.bam ${sample_id}.hisat3n_dna.multi_aligned.bam
          rm ${sample_id}.hisat3n_dna.unmapped.fastq
       }

      # define lists of r1 and r2 fq files
      R1_files=($(ls batch*/ | grep "\-R1.fq.gz"))
      R2_files=($(ls batch*/ | grep "\-R2.fq.gz"))
        
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

task hisat_single_end {
    input {
        File split_fq_tar
        File genome_fa
        File tarred_index_files
        String plate_id
        String docker

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
        mkdir -p /cromwell_root/merged_sort_bams
        mkdir -p /cromwell_root/read_overlap
   
        # define lists of r1 and r2 fq files
        R1_files=($(ls | grep "\.hisat3n_dna.split_reads.R1.fastq"))
        R2_files=($(ls | grep "\.hisat3n_dna.split_reads.R2.fastq"))

        task() {
          BASE=$(basename "$file" ".hisat3n_dna.split_reads.R1.fastq")
          echo $BASE
          echo "Running hisat on sample_id_R1" $BASE
          
          echo "Hisat 3n R1" 
          start=$(date +%s) 
   
          # hisat on R1 single end
          hisat-3n /cromwell_root/$genome_fa_basename \
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
         hisat-3n /cromwell_root/$genome_fa_basename \
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

         # remove_overlap_read_parts
         echo "call remove_overlap_read_parts" 
         start=$(date +%s) 
         python3 -c 'from cemba_data.hisat3n import *;import os;remove_overlap_read_parts(in_bam_path=os.path.join(os.path.sep,"cromwell_root","'"$BASE"'.name_sorted.filtered.bam"),out_bam_path=os.path.join(os.path.sep,"cromwell_root","'"$BASE"'.hisat3n_dna.split_reads.read_overlap.bam"))'  
         end=$(date +%s) 
         elapsed=$((end - start))  
         echo "Elapsed time to run remove overlap $elapsed seconds"
      
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
        ls 

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
  
task merge_sort_analyze {
    input {
        String plate_id
        File paired_end_unique_tar
        File read_overlap_tar
        String docker

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
      mkdir /cromwell_root/output_bams
      mkdir /cromwell_root/temp
      mkdir /cromwell_root/allc-${mcg_context}
      
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
        echo "Sort all reads by name"
        samtools sort -O BAM -@4 -m1g -o "${sample_id}.hisat3n_dna.all_reads.pos_sort.bam" "${sample_id}.hisat3n_dna.all_reads.name_sort.bam" 
        end=$(date +%s) 
        elapsed=$((end - start))  
        echo "Elapsed time to run sort by pos $elapsed seconds"
        
        start=$(date +%s)  
        echo "Call Picard remove duplicates"
        name=${sample_id}.hisat3n_dna.all_reads.deduped
        picard MarkDuplicates I=${sample_id}.hisat3n_dna.all_reads.pos_sort.bam O=/cromwell_root/output_bams/${name}.bam \
        M=/cromwell_root/output_bams/${name}.matrix.txt \
        REMOVE_DUPLICATES=true TMP_DIR=/cromwell_root/temp
        end=$(date +%s) 
        elapsed=$((end - start))  
        echo "Elapsed time to run picard $elapsed seconds"
        
        start=$(date +%s)  
        echo "Call samtools index"
        samtools index /cromwell_root/output_bams/${name}.bam
        end=$(date +%s) 
        elapsed=$((end - start)) 
        echo "Elapsed time to samtools index $elapsed seconds" 
        
        start=$(date +%s)  
        echo "Call chromatin contacts from name sorted bams" 
        python3 -c 'from cemba_data.hisat3n import *;import os;import glob;call_chromatin_contacts(bam_path="'"$sample_id"'.hisat3n_dna.all_reads.name_sort.bam",contact_prefix="'"$sample_id"'.hisat3n_dna.all_reads",save_raw=False,save_hic_format=True)'
        end=$(date +%s) 
        elapsed=$((end - start)) 
        echo "Elapsed time to chromatin contacts $elapsed seconds"

        start=$(date +%s)  
        echo "Call allcools bam-to-allc from deduped.bams" 
        /opt/conda/bin/allcools bam-to-allc \
        --bam_path /cromwell_root/output_bams/${name}.bam \
        --reference_fasta /cromwell_root/reference/~{genome_base} \
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
        --output_prefix /cromwell_root/allc-${mcg_context}/${sample_id} \
        --mc_contexts ${mcg_context} \
        --chrom_size_path ~{chromosome_sizes}
        end=$(date +%s) 
        elapsed=$((end - start)) 
        echo "Elapsed time to allcools extract-all $elapsed seconds"
        
        echo "Remove some bams"
        rm ${sample_id}.hisat3n_dna.all_reads.bam
        rm ${sample_id}.hisat3n_dna.all_reads.pos_sort.bam
        rm /cromwell_root/${sample_id}.hisat3n_dna.split_reads.read_overlap.bam
        rm /cromwell_root/${sample_id}.hisat3n_dna.unique_aligned.bam
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

      echo "Tar files."
      
      tar -cf - output_bams/*.matrix.txt | pigz > ~{plate_id}.dedup_unique_bam_and_index_unique_bam_stats.tar.gz
      tar -cf - *.hisat3n_dna.all_reads.name_sort.bam | pigz > ~{plate_id}.hisat3n_dna.all_reads.name_sort.tar.gz
    
      # tar outputs of call_chromatin_contacts
      tar -cf - *.hisat3n_dna.all_reads.3C.contact.tsv.gz | pigz > ~{plate_id}.hisat3n_dna.all_reads.3C.contact.tar.gz
      tar -cf - *.hisat3n_dna.all_reads.dedup_contacts.tsv.gz | pigz > ~{plate_id}.hisat3n_dna.all_reads.dedup_contacts.tar.gz
      tar -cf - *.hisat3n_dna.all_reads.contact_stats.csv | pigz > ~{plate_id}.chromatin_contact_stats.tar.gz
      
      # tar outputs of allcools
      tar -cf - *.allc.tsv.gz | pigz > ~{plate_id}.allc.tsv.tar.gz
      tar -cf - *.allc.tsv.gz.tbi | pigz > ~{plate_id}.allc.tbi.tar.gz
      tar -cf -  *.allc.tsv.gz.count.csv | pigz > ~{plate_id}.allc.count.tar.gz
      tar -cf -  *.tbi | pigz > ~{plate_id}.extract-allc_tbi.tar.gz
      tar -cf -  /cromwell_root/allc-${mcg_context}/*.gz | pigz > ~{plate_id}.extract-allc.tar.gz
      tar -cf -  /cromwell_root/allc-${mcg_context}/*.tbi | pigz > ~{plate_id}.extract-allc_tbi.tar.gz
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

task summary {
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

        String docker
        Int disk_size = 80
        Int mem_size = 5
        Int preemptible_tries = 3
        Int cpu = 4
    }
    command <<<
        set -euo pipefail

        mkdir /cromwell_root/fastq
        mkdir /cromwell_root/bam
        mkdir /cromwell_root/allc
        mkdir /cromwell_root/hic

        extract_and_remove() {
            if [ $# -eq 0 ];
                then
                    echo "No files exist"
                    return
            fi
            for tar in "${@}"; do
                tar -xf "$tar"
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

        mv *.trimmed.stats.txt /cromwell_root/fastq
        mv *.hisat3n_dna_summary.txt *.hisat3n_dna_split_reads_summary.R1.txt *.hisat3n_dna_split_reads_summary.R2.txt /cromwell_root/bam
        mv output_bams/*.hisat3n_dna.all_reads.deduped.matrix.txt /cromwell_root/bam
        mv *.hisat3n_dna.all_reads.contact_stats.csv /cromwell_root/hic
        mv *.allc.tsv.gz.count.csv /cromwell_root/allc
        mv cromwell_root/allc-CGN/*.allc.tsv.gz.tbi /cromwell_root/allc

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
