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
        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
        String single_end_hisat_cpu_platform = "Intel Ice Lake"
        String merge_sort_analyze_cpu_platform = "Intel Ice Lake"
    }

    # version of the pipeline
    String pipeline_version = "3.0.0"

    call Demultiplexing {
        input:
            fastq_input_read1 = fastq_input_read1,
            fastq_input_read2 = fastq_input_read2,
            random_primer_indexes = random_primer_indexes,
            plate_id = plate_id,
            batch_number = batch_number
    }

    scatter(tar in Demultiplexing.tarred_demultiplexed_fastqs) {
        call all_tasks {
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
                num_upstr_bases = num_upstr_bases,
                num_downstr_bases = num_downstr_bases,
                compress_level = compress_level,
                plate_id = plate_id,
                docker = docker
        }
    }
   
    call summary {
        input:
            trimmed_stats = all_tasks.trim_stats_tar,
            hisat3n_stats = all_tasks.hisat3n_paired_end_stats_tar,
            r1_hisat3n_stats = all_tasks.hisat3n_dna_split_reads_summary_R1_tar,
            r2_hisat3n_stats = all_tasks.hisat3n_dna_split_reads_summary_R2_tar,
            plate_id = plate_id
    }


    output {
        File MappingSummary = summary.mapping_summary
        Array[File] name_sorted_bams = all_tasks.name_sorted_bam
        Array[File] unique_reads_cgn_extraction_allc= all_tasks.allc
        Array[File] unique_reads_cgn_extraction_tbi = all_tasks.tbi
        Array[File] reference_version = all_tasks.reference_version
        Array[File] all_reads_dedup_contacts = all_tasks.all_reads_dedup_contacts
        Array[File] all_reads_3C_contacts = all_tasks.all_reads_3C_contacts
        Array[File] chromatin_contact_stats = all_tasks.chromatin_contact_stats
        Array[File] unique_reads_cgn_extraction_allc_extract = all_tasks.extract_allc_output_allc_tar
        Array[File] unique_reads_cgn_extraction_tbi_extract = all_tasks.extract_allc_output_tbi_tar

    }
}

task Demultiplexing {
  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_id
    Int batch_number

    String docker_image = "us.gcr.io/broad-gotc-prod/hisat3n:2.1.0-2.2.1-1709740155"
    Int disk_size = 1000
    Int mem_size = 10
    Int preemptible_tries = 3
    Int cpu = 8
  }

  command <<<
    set -euo pipefail

    # Cat files for each r1, r2
    cat ~{sep=' ' fastq_input_read1} > r1.fastq.gz
    cat ~{sep=' ' fastq_input_read2} > r2.fastq.gz

    /opt/conda/bin/cutadapt -Z -e 0.01 --no-indels \
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
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: cpu
    memory: "${mem_size} GiB"
    preemptible: preemptible_tries
  }

  output {
    Array[File] tarred_demultiplexed_fastqs = glob("*.tar.gz")
    File stats = "~{plate_id}.stats.txt"
    }
}


task all_tasks{
    input {
        File tarred_demultiplexed_fastqs
        File tarred_index_files
       
        File genome_fa
        String genome_base = basename(genome_fa)
        File chromosome_sizes
        String plate_id

        # hisat paired task -- cutadapt 
        String r1_adapter
        String r2_adapter
        Int r1_left_cut
        Int r1_right_cut
        Int r2_left_cut
        Int r2_right_cut
        Int min_read_length

        # merge, sort and analyze
        Int num_upstr_bases
        Int num_downstr_bases
        Int compress_level

        # run time variables
        Int disk_size = 2000
        Int cpu = 64
        Int mem_size = 128
        String cpu_platform = "Intel Ice Lake"
        String docker
    }

    command <<<
        set -euo pipefail
        set -x
        lscpu
  
        min_read_length=~{min_read_length}

        # check genomic reference version and print to output txt file
        STRING=~{genome_fa}
        BASE=$(basename $STRING .fa)

        # save reference version in seperate file 
        echo "The reference is $BASE" > ~{plate_id}.reference_version.txt

        # untar the index files for hisat task
        start=$(date +%s)
        echo "Untarring tarred_index_files"
        tar -zxvf ~{tarred_index_files}
        rm ~{tarred_index_files}
        end=$(date +%s) 
        elapsed=$((end - start)) 
        echo "Elapsed time to untar tarred_index_files: $elapsed seconds"

        # reference and index 
        start=$(date +%s)  
        echo "Reference and index fasta"
        cp ~{genome_fa} .
        samtools faidx $BASE.fa
        ls
        end=$(date +%s) 
        elapsed=$((end - start)) 
        echo "Elapsed time to index fasta $elapsed seconds"
 
        # untar the demultiplexed fastqs for sort and trim task
        start=$(date +%s)  
        echo "Reference and index fasta"
        tar -xf ~{tarred_demultiplexed_fastqs}
        end=$(date +%s) 
        elapsed=$((end - start)) 
        echo "Elapsed time to index fasta $elapsed seconds"
      

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
        
        paired_end_hisat() {
          local sample_id=$1
          echo $sample_id
          
          r2_file="${sample_id}-R2.fq.gz"
          r1_file="${sample_id}-R1.fq.gz"
          
          # sort R1 
          start=$(date +%s)
          echo "Run sort r1"
          zcat /cromwell_root/batch*/"$r1_file" | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > "${sample_id}-R1_sorted.fq"
          end=$(date +%s) 
          elapsed=$((end - start)) 
          echo "Elapsed time to run sort r1: $elapsed seconds"
    
          # sort R2
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
          --report=minimal \
          -O 6 \
          -q 20 \
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
          hisat-3n /cromwell_root/$BASE \
          -q \
          -1 ${sample_id}-R1_trimmed.fq.gz \
          -2 ${sample_id}-R2_trimmed.fq.gz \
          --directional-mapping-reverse \
          --base-change C,T \
          --repeat \
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
    
          # why? wait 15 seconds for the files to be written
          # sleep 15

          # remove files
          rm /cromwell_root/batch*/${sample_id}-R1.fq.gz /cromwell_root/batch*/${sample_id}-R2.fq.gz
          rm ${sample_id}-R1_sorted.fq ${sample_id}-R2_sorted.fq
          rm ${sample_id}-R1_trimmed.fq.gz ${sample_id}-R2_trimmed.fq.gz
          rm ${sample_id}.hisat3n_dna.unsort.bam ${sample_id}.hisat3n_dna.multi_aligned.bam
          rm ${sample_id}.hisat3n_dna.unmapped.fastq

        }        

        single_end_hisat() {
          local sample_id=$1
          echo $sample_id
          echo "Running hisat on sample_id_R1" $sample_id
          
          echo "Hisat single end R1" 
          start=$(date +%s) 
   
          # hisat on R1 single end hisat3n_dna.split_reads.R1.fastq
          hisat-3n /cromwell_root/$BASE \
          -q \
          -U ${sample_id}.hisat3n_dna.split_reads.R1.fastq \
          -S ${sample_id}.hisat3n_dna.split_reads.R1.sam --directional-mapping-reverse --base-change C,T \
          --repeat \
          --no-spliced-alignment \
          --no-temp-splicesite \
          -t \
          --new-summary \
          --summary-file ${sample_id}.hisat3n_dna_split_reads_summary.R1.txt \
          --threads 8 
        
         end=$(date +%s) 
         elapsed=$((end - start))  
         echo "Elapsed time to run $elapsed seconds"
         echo "Finish running hisat on sample_id_R1" $sample_id
         
         echo "Hisat single end R2" 
         start=$(date +%s)        
         echo "Running hisat on sample_id_R2" $sample_id

         # hisat on R2 single end
         hisat-3n /cromwell_root/$BASE \
         -q \
         -U ${sample_id}.hisat3n_dna.split_reads.R2.fastq \
         -S ${sample_id}.hisat3n_dna.split_reads.R2.sam --directional-mapping --base-change C,T \
         --repeat \
         --no-spliced-alignment \
         --no-temp-splicesite \
         -t --new-summary \
         --summary-file ${sample_id}.hisat3n_dna_split_reads_summary.R2.txt \
         --threads 8

         end=$(date +%s) 
         elapsed=$((end - start)) 
         echo "Elapsed time to run $elapsed seconds"
         echo "Finish running hisat on sample_id_R2" $sample_id
        
         # samtools merge
         echo "samtools merge R1 and R2" 
         start=$(date +%s)        
         samtools merge -o ${sample_id}.name_merged.sam ${sample_id}.hisat3n_dna.split_reads.R1.sam ${sample_id}.hisat3n_dna.split_reads.R2.sam -@8
         end=$(date +%s) 
         elapsed=$((end - start))  
         echo "Elapsed time to run samtools merge $elapsed seconds"
                  
         # samtools sort 
         echo "samtools sort R1 and R2" 
         start=$(date +%s)        
         samtools sort -n -@8 -m1g ${sample_id}.name_merged.sam -o ${sample_id}.name_sorted.bam
         end=$(date +%s) 
         elapsed=$((end - start)) 
         echo "Elapsed time to run samtools sort $elapsed seconds"

         # samtools filter bam
         echo "samtools -q 10" 
         start=$(date +%s)  
         samtools view -q 10 ${sample_id}.name_sorted.bam -o ${sample_id}.name_sorted.filtered.bam
         end=$(date +%s) 
         elapsed=$((end - start)) 
         echo "Elapsed time to run samtools -q 10 $elapsed seconds"

         # remove_overlap_read_parts
         echo "call remove_overlap_read_parts" 
         start=$(date +%s) 
         python3 -c 'from cemba_data.hisat3n import *;import os;remove_overlap_read_parts(in_bam_path=os.path.join(os.path.sep,"cromwell_root","'"$sample_id"'.name_sorted.filtered.bam"),out_bam_path=os.path.join(os.path.sep,"cromwell_root","'"$sample_id"'.hisat3n_dna.split_reads.read_overlap.bam"))'  
         end=$(date +%s) 
         elapsed=$((end - start))  
         echo "Elapsed time to run remove overlap $elapsed seconds"

         # remove files? what files?
         rm ${sample_id}.name_sorted.bam
         rm ${sample_id}.name_merged.sam
         rm ${sample_id}.hisat3n_dna.split_reads.R1.sam
         rm ${sample_id}.hisat3n_dna.split_reads.R2.sam
         rm ${sample_id}.name_sorted.filtered.bam #needed? 
        }

        merge_sort_analyze() {
          local file=$1
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
          --reference_fasta /cromwell_root/~{genome_base} \
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
        
        # define lists of r1 and r2 fq files
        R1_files=($(ls batch*/ | grep "\-R1.fq.gz"))
        R2_files=($(ls batch*/ | grep "\-R2.fq.gz"))

        # run 4 instances in parallel each with 8 threads
        for file in "${R1_files[@]}"; do
          (
            echo "starting $file.."
            # nova_seq_full-A1-R1.fq.gz
            sample_id=$(basename "$file" "-R1.fq.gz")
            echo $sample_id
            paired_end_hisat "$sample_id"
            single_end_hisat "$sample_id"
            merge_sort_analyze "$sample_id"
            sleep $(( (RANDOM % 3) + 1))
          ) &
          if [[ $(jobs -r -p | wc -l) -ge 8 ]]; then
            wait -n
          fi
        done

        wait
        echo "All done running tasks."
        
        echo "Tar up files"
        echo "Tar up summary files for trimmming + paired end"
        start=$(date +%s)
        tar -cf - *.trimmed.stats.txt | pigz > ~{plate_id}.trimmed_stats_files.tar.gz
        tar -cf - *.hisat3n_dna_summary.txt | pigz > ~{plate_id}.hisat3n_paired_end_stats_files.tar.gz
        end=$(date +%s) 
        elapsed=$((end - start))  
        echo "Elapsed time to run tar summary files for trimming + paired end $elapsed seconds"

        echo "Tar up summary files for single end"
        start=$(date +%s)
        tar -cf - *.hisat3n_dna_split_reads_summary.R1.txt | pigz > ~{plate_id}.hisat3n_dna_split_reads_summary.R1.tar.gz
        tar -cf - *.hisat3n_dna_split_reads_summary.R2.txt | pigz > ~{plate_id}.hisat3n_dna_split_reads_summary.R2.tar.gz
        end=$(date +%s) 
        elapsed=$((end - start))  
        echo "Elapsed time to run tar summary files for single end $elapsed seconds"
        
        echo "Tar up bams (*all_reads.name_sort.bam) files"
        start=$(date +%s)
        tar -cf - *.hisat3n_dna.all_reads.name_sort.bam | pigz > ~{plate_id}.hisat3n_dna.all_reads.name_sort.tar.gz
        end=$(date +%s) 
        elapsed=$((end - start))  
        echo "Elapsed time to run tar name_sort.bam files $elapsed seconds"
        
        # tar outputs of call_chromatin_contacts
        echo "Tar up 3C.contact.tsv.gz, dedup_contacts.tsv.gz and contact_stats.csv"
        start=$(date +%s)
        tar -cf - *.hisat3n_dna.all_reads.3C.contact.tsv.gz | pigz > ~{plate_id}.hisat3n_dna.all_reads.3C.contact.tar.gz
        tar -cf - *.hisat3n_dna.all_reads.dedup_contacts.tsv.gz | pigz > ~{plate_id}.hisat3n_dna.all_reads.dedup_contacts.tar.gz
        tar -cf - *.hisat3n_dna.all_reads.contact_stats.csv | pigz > ~{plate_id}.chromatin_contact_stats.tar.gz
        end=$(date +%s) 
        elapsed=$((end - start))  
        echo "Elapsed time to run tar name_sort.bam files $elapsed seconds"
        
        # tar outputs of allcools
        echo "Tar up allc.tsv, .allc.tbi, allc.count... "
        start=$(date +%s)  
        tar -cf - *.allc.tsv.gz | pigz > ~{plate_id}.allc.tsv.tar.gz
        tar -cf - *.allc.tsv.gz.tbi | pigz > ~{plate_id}.allc.tbi.tar.gz
        tar -cf - *.allc.tsv.gz.count.csv | pigz > ~{plate_id}.allc.count.tar.gz
        #tar -cf - *.tbi | pigz > ~{plate_id}.extract-allc_tbi.tar.gz ? -- what is this -- same name as before
        tar -cf - /cromwell_root/allc-${mcg_context}/*.gz | pigz > ~{plate_id}.extract-allc.tar.gz
        tar -cf - /cromwell_root/allc-${mcg_context}/*.tbi | pigz > ~{plate_id}.extract-allc_tbi.tar.gz
        end=$(date +%s) 
        elapsed=$((end - start))  
        echo "Elapsed time to run tar files $elapsed seconds"    

        du -h *
    >>>

    runtime {
        docker: docker
        disks: "local-disk ${disk_size} SSD"
        cpu: cpu
        memory: "${mem_size} GiB"
        cpuPlatform: cpu_platform
    }

    output {
        # bams
        File name_sorted_bam = "~{plate_id}.hisat3n_dna.all_reads.name_sort.tar.gz"
        # summary outputs for paired/single end
        File trim_stats_tar = "~{plate_id}.trimmed_stats_files.tar.gz"
        File hisat3n_paired_end_stats_tar = "~{plate_id}.hisat3n_paired_end_stats_files.tar.gz"
        File hisat3n_dna_split_reads_summary_R1_tar = "~{plate_id}.hisat3n_dna_split_reads_summary.R1.tar.gz"
        File hisat3n_dna_split_reads_summary_R2_tar = "~{plate_id}.hisat3n_dna_split_reads_summary.R2.tar.gz"
        # chromatin contacts
        File all_reads_dedup_contacts = "~{plate_id}.hisat3n_dna.all_reads.dedup_contacts.tar.gz"
        File all_reads_3C_contacts = "~{plate_id}.hisat3n_dna.all_reads.3C.contact.tar.gz"
        File chromatin_contact_stats = "~{plate_id}.chromatin_contact_stats.tar.gz"
        # allc 
        File allc = "~{plate_id}.allc.tsv.tar.gz"
        File tbi = "~{plate_id}.allc.tbi.tar.gz"
        File allc_uniq_reads_stats = "~{plate_id}.allc.count.tar.gz"
        File extract_allc_output_tbi_tar = "~{plate_id}.extract-allc_tbi.tar.gz"
        File extract_allc_output_allc_tar  = "~{plate_id}.extract-allc.tar.gz"
        # reference version
        File reference_version = "~{plate_id}.reference_version.txt"    
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

        String docker = "us.gcr.io/broad-gotc-prod/hisat3n:2.1.0-2.2.1-1709740155"
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