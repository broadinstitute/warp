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

    }

    # version of the pipeline
    String pipeline_version = "2.0.1"

    call Demultiplexing {
        input:
            fastq_input_read1 = fastq_input_read1,
            fastq_input_read2 = fastq_input_read2,
            random_primer_indexes = random_primer_indexes,
            plate_id = plate_id,
            batch_number = batch_number
    }

    scatter(tar in Demultiplexing.tarred_demultiplexed_fastqs) {
        call Sort_and_trim_r1_and_r2 {
            input:
                tarred_demultiplexed_fastqs = tar,
                r1_adapter = r1_adapter,
                r2_adapter = r2_adapter,
                r1_left_cut = r1_left_cut,
                r1_right_cut = r1_right_cut,
                r2_left_cut = r2_left_cut,
                r2_right_cut = r2_right_cut,
                min_read_length = min_read_length,
                plate_id = plate_id
        }

        call Hisat_3n_pair_end_mapping_dna_mode {
            input:
                r1_trimmed_tar = Sort_and_trim_r1_and_r2.r1_trimmed_fq_tar,
                r2_trimmed_tar = Sort_and_trim_r1_and_r2.r2_trimmed_fq_tar,
                tarred_index_files = tarred_index_files,
                genome_fa = genome_fa,
                chromosome_sizes = chromosome_sizes,
                plate_id = plate_id
        }

        call Separate_and_split_unmapped_reads {
            input:
                hisat3n_bam_tar = Hisat_3n_pair_end_mapping_dna_mode.hisat3n_paired_end_bam_tar,
                min_read_length = min_read_length,
                plate_id = plate_id,
        }

        call Hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name_and_remove_overlap {
            input:
                split_fq_tar = Separate_and_split_unmapped_reads.split_fq_tar,
                tarred_index_files = tarred_index_files,
                genome_fa = genome_fa,
                plate_id = plate_id
        }

        call merge_original_and_split_bam_and_sort_all_reads_by_name_and_position_and_deduplicate {
            input:
                bam = Separate_and_split_unmapped_reads.unique_bam_tar,
                split_bam = Hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name_and_remove_overlap.remove_overlaps_output_bam_tar,
                plate_id = plate_id
        }

        call call_chromatin_contacts {
            input:
                name_sorted_bam = merge_original_and_split_bam_and_sort_all_reads_by_name_and_position_and_deduplicate.name_sorted_bam,
                plate_id = plate_id
        }

        call unique_reads_allc_and_cgn_extraction {
            input:
                bam_and_index_tar = merge_original_and_split_bam_and_sort_all_reads_by_name_and_position_and_deduplicate.dedup_output_bam_tar,
                genome_fa = genome_fa,
                num_upstr_bases = num_upstr_bases,
                num_downstr_bases = num_downstr_bases,
                compress_level = compress_level,
                plate_id = plate_id,
                chromosome_sizes = chromosome_sizes
        }
    }

    call summary {
        input:
            trimmed_stats = Sort_and_trim_r1_and_r2.trim_stats_tar,
            hisat3n_stats = Hisat_3n_pair_end_mapping_dna_mode.hisat3n_paired_end_stats_tar,
            r1_hisat3n_stats = Hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name_and_remove_overlap.hisat3n_dna_split_reads_summary_R1_tar,
            r2_hisat3n_stats = Hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name_and_remove_overlap.hisat3n_dna_split_reads_summary_R2_tar,
            dedup_stats = merge_original_and_split_bam_and_sort_all_reads_by_name_and_position_and_deduplicate.dedup_stats_tar,
            chromatin_contact_stats = call_chromatin_contacts.chromatin_contact_stats,
            allc_uniq_reads_stats = unique_reads_allc_and_cgn_extraction.allc_uniq_reads_stats,
            unique_reads_cgn_extraction_tbi = unique_reads_allc_and_cgn_extraction.extract_allc_output_tbi_tar,
            plate_id = plate_id
    }

    output {
        File MappingSummary = summary.mapping_summary
        Array[File] name_sorted_bams = merge_original_and_split_bam_and_sort_all_reads_by_name_and_position_and_deduplicate.name_sorted_bam
        Array[File] unique_reads_cgn_extraction_allc= unique_reads_allc_and_cgn_extraction.allc
        Array[File] unique_reads_cgn_extraction_tbi = unique_reads_allc_and_cgn_extraction.tbi
        Array[File] unique_reads_cgn_extraction_allc_extract = unique_reads_allc_and_cgn_extraction.extract_allc_output_allc_tar
        Array[File] unique_reads_cgn_extraction_tbi_extract = unique_reads_allc_and_cgn_extraction.extract_allc_output_tbi_tar
        Array[File] reference_version = Hisat_3n_pair_end_mapping_dna_mode.reference_version
        Array[File] chromatin_contact_stats = call_chromatin_contacts.chromatin_contact_stats
        Array[File] all_reads_dedup_contacts = call_chromatin_contacts.all_reads_dedup_contacts
        Array[File] all_reads_3C_contacts = call_chromatin_contacts.all_reads_3C_contacts
    }
}

task Demultiplexing {
  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_id
    Int batch_number

    String docker_image = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
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

task Sort_and_trim_r1_and_r2 {
    input {
        File tarred_demultiplexed_fastqs
        String plate_id
        String r1_adapter
        String r2_adapter
        Int r1_left_cut
        Int r1_right_cut
        Int r2_left_cut
        Int r2_right_cut
        Int min_read_length

        Int disk_size = 500
        Int mem_size = 16
        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
        Int preemptible_tries = 3
        Int cpu = 4

    }
    command <<<
    set -euo pipefail

    # untar the demultiplexed fastqs
    tar -xf ~{tarred_demultiplexed_fastqs}

    #change into batch subfolder
    cd batch*
    # define lists of r1 and r2 fq files
    R1_files=($(ls | grep "\-R1.fq.gz"))
    R2_files=($(ls | grep "\-R2.fq.gz"))

    # loop over R1 and R2 files and sort them
    for file in "${R1_files[@]}"; do
      sample_id=$(basename "$file" "-R1.fq.gz")
      r2_file="${sample_id}-R2.fq.gz"
      zcat "$file" | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > "${sample_id}-R1_sorted.fq"
      zcat "$r2_file" | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > "${sample_id}-R2_sorted.fq"
    done


    echo "Starting to trim with Cutadapt"
    sorted_R1_files=($(ls | grep "\-R1_sorted.fq"))
    for file in "${sorted_R1_files[@]}"; do
      sample_id=$(basename "$file" "-R1_sorted.fq")
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
        -m ~{min_read_length}:~{min_read_length} \
        --pair-filter 'both' \
        -o ${sample_id}-R1_trimmed.fq.gz \
        -p ${sample_id}-R2_trimmed.fq.gz \
        ${sample_id}-R1_sorted.fq ${sample_id}-R2_sorted.fq \
        > ${sample_id}.trimmed.stats.txt
    done

    echo "Tarring up the trimmed files and stats files"

    tar -zcvf ~{plate_id}.R1_trimmed_files.tar.gz *-R1_trimmed.fq.gz
    tar -zcvf ~{plate_id}.R2_trimmed_files.tar.gz *-R2_trimmed.fq.gz
    tar -zcvf ~{plate_id}.trimmed_stats_files.tar.gz *.trimmed.stats.txt
    # move files back to root
    mv ~{plate_id}.R1_trimmed_files.tar.gz ../~{plate_id}.R1_trimmed_files.tar.gz
    mv ~{plate_id}.R2_trimmed_files.tar.gz ../~{plate_id}.R2_trimmed_files.tar.gz
    mv ~{plate_id}.trimmed_stats_files.tar.gz ../~{plate_id}.trimmed_stats_files.tar.gz
    >>>
    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
    }
    output {
        File r1_trimmed_fq_tar = "~{plate_id}.R1_trimmed_files.tar.gz"
        File r2_trimmed_fq_tar = "~{plate_id}.R2_trimmed_files.tar.gz"
        File trim_stats_tar = "~{plate_id}.trimmed_stats_files.tar.gz"
    }
}

task Hisat_3n_pair_end_mapping_dna_mode{
    input {
        File r1_trimmed_tar
        File r2_trimmed_tar
        File tarred_index_files
        File genome_fa
        File chromosome_sizes
        String plate_id

        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
        Int disk_size = 1000
        Int mem_size = 64
        Int preemptible_tries = 3
        Int cpu = 16
    }
    command <<<
        set -euo pipefail

        # check genomic reference version and print to output txt file
        STRING=~{genome_fa}
        BASE=$(basename $STRING .fa)

        echo "The reference is $BASE" > ~{plate_id}.reference_version.txt

        # untar the index files
        echo "Untarring the index files"
        tar -zxvf ~{tarred_index_files}
        rm ~{tarred_index_files}

        cp ~{genome_fa} .

        #get the basename of the genome_fa file
        genome_fa_basename=$(basename ~{genome_fa} .fa)
        echo "samtools faidx $genome_fa_basename.fa"
        samtools faidx $genome_fa_basename.fa

        # untar the demultiplexed fastq files
        echo "Untarring the fastq files"
        tar -zxvf ~{r1_trimmed_tar}
        tar -zxvf ~{r2_trimmed_tar}
        rm ~{r1_trimmed_tar}
        rm ~{r2_trimmed_tar}

        # define lists of r1 and r2 fq files
        R1_files=($(ls | grep "\-R1_trimmed.fq.gz"))
        R2_files=($(ls | grep "\-R2_trimmed.fq.gz"))

        for file in "${R1_files[@]}"; do
          sample_id=$(basename "$file" "-R1_trimmed.fq.gz")
          hisat-3n /cromwell_root/$genome_fa_basename \
          -q \
          -1 ${sample_id}-R1_trimmed.fq.gz \
          -2 ${sample_id}-R2_trimmed.fq.gz \
          --directional-mapping-reverse \
          --base-change C,T \
          --no-repeat-index \
          --no-spliced-alignment \
          --no-temp-splicesite \
          -t \
          --new-summary \
          --summary-file ${sample_id}.hisat3n_dna_summary.txt \
          --threads 11 | samtools view -b -q 0 -o "${sample_id}.hisat3n_dna.unsort.bam"
        done

        # tar up the bam files and stats files
        tar -zcvf ~{plate_id}.hisat3n_paired_end_bam_files.tar.gz *.bam
        tar -zcvf ~{plate_id}.hisat3n_paired_end_stats_files.tar.gz *.hisat3n_dna_summary.txt


    >>>
    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
    }
    output {
        File hisat3n_paired_end_bam_tar = "~{plate_id}.hisat3n_paired_end_bam_files.tar.gz"
        File hisat3n_paired_end_stats_tar = "~{plate_id}.hisat3n_paired_end_stats_files.tar.gz"
        File reference_version = "~{plate_id}.reference_version.txt"
    }
}

task Separate_and_split_unmapped_reads {
    input {
        File hisat3n_bam_tar
        Int min_read_length
        String plate_id

        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
        Int disk_size = 1000
        Int mem_size = 10
        Int preemptible_tries = 3
        Int cpu = 8

    }
    command <<<

        set -euo pipefail

        # untar the hisat3n bam files
        tar -xf ~{hisat3n_bam_tar}
        rm ~{hisat3n_bam_tar}

        python3 <<CODE

        from cemba_data.hisat3n import separate_unique_and_multi_align_reads
        import os
        import glob

        pattern = "*.hisat3n_dna.unsort.bam"
        bam_files = glob.glob(os.path.join('/cromwell_root/', pattern))


        for file in bam_files:
            full_filename = os.path.basename(file)
            sample_id = full_filename.replace(".hisat3n_dna.unsort.bam", "")
            in_bam_path = f"{sample_id}.hisat3n_dna.unsort.bam"
            out_unique_path = f"{sample_id}.hisat3n_dna.unique_aligned.bam"
            out_multi_path = f"{sample_id}.hisat3n_dna.multi_aligned.bam"
            out_unmappable_path = f"{sample_id}.hisat3n_dna.unmapped.fastq"

            separate_unique_and_multi_align_reads(
                in_bam_path=in_bam_path,
                out_unique_path=out_unique_path,
                out_multi_path=out_multi_path,
                out_unmappable_path=out_unmappable_path,
                unmappable_format='fastq',
                mapq_cutoff=10,
                qlen_cutoff=~{min_read_length}
            )

        CODE

        # tar up the uniqe bams
        tar -zcvf ~{plate_id}.hisat3n_paired_end_unique_bam_files.tar.gz *.hisat3n_dna.unique_aligned.bam

        # tar up the unmapped fastq files
        tar -zcvf ~{plate_id}.hisat3n_paired_end_unmapped_fastq_files.tar.gz *.hisat3n_dna.unmapped.fastq

        # untar the unmapped fastq files
        tar -xf ~{plate_id}.hisat3n_paired_end_unmapped_fastq_files.tar.gz

        python3 <<CODE

        from cemba_data.hisat3n import *
        import os
        import glob

        pattern = "*.hisat3n_dna.unmapped.fastq"
        fastq_files = glob.glob(os.path.join('/cromwell_root/', pattern))

        for file in fastq_files:
          full_filename = os.path.basename(file)
          sample_id = full_filename.replace(".hisat3n_dna.unmapped.fastq", "")
          fastq_path = f"{sample_id}.hisat3n_dna.unmapped.fastq"
          output_prefix = f"{sample_id}.hisat3n_dna.split_reads"

          split_hisat3n_unmapped_reads(
            fastq_path=fastq_path,
            output_prefix=output_prefix,
            min_length=~{min_read_length}
          )

        CODE

        # wait 15 seconds for the files to be written
        sleep 15

        # tar up the split fastq files
        tar -zcvf ~{plate_id}.hisat3n_paired_end_split_fastq_files.tar.gz *.split_reads*.fastq

    >>>
    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
    }
    output {
        File unique_bam_tar = "~{plate_id}.hisat3n_paired_end_unique_bam_files.tar.gz"
        File split_fq_tar = "~{plate_id}.hisat3n_paired_end_split_fastq_files.tar.gz"
    }
}

task Hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name_and_remove_overlap {
    input {
        File split_fq_tar
        File genome_fa
        File tarred_index_files
        String plate_id

        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
        Int disk_size = 500
        Int mem_size = 64
        Int preemptible_tries = 3
        Int cpu = 16
    }
    command <<<
        set -euo pipefail


        # untar the tarred index files
        tar -xvf ~{tarred_index_files}
        rm ~{tarred_index_files}

        cp ~{genome_fa} .

        #get the basename of the genome_fa file
        genome_fa_basename=$(basename ~{genome_fa} .fa)
        samtools faidx $genome_fa_basename.fa

        # untar the unmapped fastq files
        tar -xvf ~{split_fq_tar}
        rm ~{split_fq_tar}

        # define lists of r1 and r2 fq files
        R1_files=($(ls | grep "\.hisat3n_dna.split_reads.R1.fastq"))
        R2_files=($(ls | grep "\.hisat3n_dna.split_reads.R2.fastq"))

        for file in "${R1_files[@]}"; do
          sample_id=$(basename "$file" ".hisat3n_dna.split_reads.R1.fastq")
          hisat-3n /cromwell_root/$genome_fa_basename \
          -q \
          -U ${sample_id}.hisat3n_dna.split_reads.R1.fastq \
          --directional-mapping-reverse \
          --base-change C,T \
          --no-repeat-index \
          --no-spliced-alignment \
          --no-temp-splicesite \
          -t \
          --new-summary \
          --summary-file ${sample_id}.hisat3n_dna_split_reads_summary.R1.txt \
          --threads 11 | samtools view -b -q 10 -o "${sample_id}.hisat3n_dna.split_reads.R1.bam"
        done

       for file in "${R2_files[@]}"; do
         sample_id=$(basename "$file" ".hisat3n_dna.split_reads.R2.fastq")
         hisat-3n /cromwell_root/$genome_fa_basename \
         -q \
         -U ${sample_id}.hisat3n_dna.split_reads.R2.fastq \
         --directional-mapping \
         --base-change C,T \
         --no-repeat-index \
         --no-spliced-alignment \
         --no-temp-splicesite \
         -t --new-summary \
         --summary-file ${sample_id}.hisat3n_dna_split_reads_summary.R2.txt \
         --threads 11 | samtools view -b -q 10 -o "${sample_id}.hisat3n_dna.split_reads.R2.bam"
       done

       # tar up the r1 and r2 stats files
       tar -zcvf ~{plate_id}.hisat3n_dna_split_reads_summary.R1.tar.gz *.hisat3n_dna_split_reads_summary.R1.txt
       tar -zcvf ~{plate_id}.hisat3n_dna_split_reads_summary.R2.tar.gz *.hisat3n_dna_split_reads_summary.R2.txt


       # define lists of r1 and r2 bam files
       R1_bams=($(ls | grep "\.hisat3n_dna.split_reads.R1.bam"))
       R2_bams=($(ls | grep "\.hisat3n_dna.split_reads.R2.bam"))

       # Loop through the R1 BAM files
       for r1_bam in "${R1_bams[@]}"; do
         # Extract the corresponding R2 BAM file
         r2_bam="${r1_bam/.hisat3n_dna.split_reads.R1.bam/.hisat3n_dna.split_reads.R2.bam}"

         # Define the output BAM file name
         output_bam="$(basename ${r1_bam/.hisat3n_dna.split_reads.R1.bam/.hisat3n_dna.split_reads.name_sort.bam})"

         # Perform the samtools merge and sort commands
         samtools merge -o - "$r1_bam" "$r2_bam" | samtools sort -n -o "$output_bam" -
       done

       #tar up the merged bam files
       tar -zcvf ~{plate_id}.hisat3n_dna.split_reads.name_sort.bam.tar.gz *.hisat3n_dna.split_reads.name_sort.bam

       # unzip bam file
       tar -xf  ~{plate_id}.hisat3n_dna.split_reads.name_sort.bam.tar.gz

       # create output dir
       mkdir /cromwell_root/output_bams
       # get bams
       bams=($(ls | grep "sort.bam$"))

       # loop through bams and run python script on each bam
       # scatter instead of for loop to optimize
       python3 <<CODE
       from cemba_data.hisat3n import *
       import os
       bams="${bams[@]}"
       for bam in bams.split(" "):
           name=".".join(bam.split(".")[:3])+".read_overlap.bam"
           remove_overlap_read_parts(in_bam_path=os.path.join(os.path.sep, "cromwell_root", bam), out_bam_path=os.path.join(os.path.sep, "cromwell_root", "output_bams", name))
       CODE

       cd /cromwell_root/output_bams

       #tar up the merged bam files
       tar -zcvf ../~{plate_id}.remove_overlap_read_parts.tar.gz *bam

    >>>
    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
    }
    output {
        #File merge_sorted_bam_tar = "~{plate_id}.hisat3n_dna.split_reads.name_sort.bam.tar.gz"
        File hisat3n_dna_split_reads_summary_R1_tar = "~{plate_id}.hisat3n_dna_split_reads_summary.R1.tar.gz"
        File hisat3n_dna_split_reads_summary_R2_tar = "~{plate_id}.hisat3n_dna_split_reads_summary.R2.tar.gz"
        File remove_overlaps_output_bam_tar = "~{plate_id}.remove_overlap_read_parts.tar.gz"
    }
}

task merge_original_and_split_bam_and_sort_all_reads_by_name_and_position_and_deduplicate {
    input {
        File bam
        File split_bam
        String plate_id

        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
        Int disk_size = 1000
        Int mem_size = 50
        Int preemptible_tries = 3
        Int cpu = 8
    }
    command <<<
      set -euo pipefail
      #unzip bam file
      tar -xf ~{bam}
      tar -xf ~{split_bam}
      rm ~{bam}
      rm ~{split_bam}

      echo "samtools merge and sort"
      # define lists of r1 and r2 fq files
      UNIQUE_BAMS=($(ls | grep "\.hisat3n_dna.unique_aligned.bam"))
      SPLIT_BAMS=($(ls | grep "\.hisat3n_dna.split_reads.read_overlap.bam"))

      for file in "${UNIQUE_BAMS[@]}"; do
        sample_id=$(basename "$file" ".hisat3n_dna.unique_aligned.bam")
        samtools merge -f "${sample_id}.hisat3n_dna.all_reads.bam" "${sample_id}.hisat3n_dna.unique_aligned.bam" "${sample_id}.hisat3n_dna.split_reads.read_overlap.bam"
        samtools sort -n -o "${sample_id}.hisat3n_dna.all_reads.name_sort.bam" "${sample_id}.hisat3n_dna.all_reads.bam"
        samtools sort -O BAM -o "${sample_id}.hisat3n_dna.all_reads.pos_sort.bam" "${sample_id}.hisat3n_dna.all_reads.name_sort.bam"
      done

      echo "Zip files"
      #tar up the merged bam files
      tar -zcvf ~{plate_id}.hisat3n_dna.all_reads.pos_sort.tar.gz *.hisat3n_dna.all_reads.pos_sort.bam
      tar -zcvf ~{plate_id}.hisat3n_dna.all_reads.name_sort.tar.gz *.hisat3n_dna.all_reads.name_sort.bam


      # unzip files
      tar -xf ~{plate_id}.hisat3n_dna.all_reads.pos_sort.tar.gz

      # create output dir
      mkdir /cromwell_root/output_bams
      mkdir /cromwell_root/temp

      # name : AD3C_BA17_2027_P1-1-B11-G13.hisat3n_dna.all_reads.pos_sort.bam
      for file in *.pos_sort.bam
      do
        name=`echo $file | cut -d. -f1`
        name=$name.hisat3n_dna.all_reads.deduped
        echo $name
        echo "Call Picard"
        picard MarkDuplicates I=$file O=/cromwell_root/output_bams/$name.bam \
        M=/cromwell_root/output_bams/$name.matrix.txt \
        REMOVE_DUPLICATES=true TMP_DIR=/cromwell_root/temp
        echo "Call samtools index"
        samtools index /cromwell_root/output_bams/$name.bam
      done

      cd /cromwell_root

      #tar up the output files
      tar -zcvf ~{plate_id}.dedup_unique_bam_and_index_unique_bam.tar.gz output_bams

      #tar up the stats files
      tar -zcvf ~{plate_id}.dedup_unique_bam_and_index_unique_bam_stats.tar.gz output_bams/*.matrix.txt

    >>>
    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
    }
    output {
        File name_sorted_bam = "~{plate_id}.hisat3n_dna.all_reads.name_sort.tar.gz"
        File dedup_output_bam_tar = "~{plate_id}.dedup_unique_bam_and_index_unique_bam.tar.gz"
        File dedup_stats_tar = "~{plate_id}.dedup_unique_bam_and_index_unique_bam_stats.tar.gz"
    }
}

task call_chromatin_contacts {
    input {
        File name_sorted_bam
        String plate_id

        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
        Int disk_size = 500
        Int mem_size = 32
        Int preemptible_tries = 3
        Int cpu = 8
    }
    command <<<
        set -euo pipefail

        # untar the name sorted bam files
        tar -xf ~{name_sorted_bam}
        rm ~{name_sorted_bam}

        python3 <<CODE

        from cemba_data.hisat3n import *
        import os
        import glob

        pattern = "*.hisat3n_dna.all_reads.name_sort.bam"
        bam_files = glob.glob(os.path.join('/cromwell_root/', pattern))

        for file in bam_files:
            full_filename = os.path.basename(file)
            sample_id = full_filename.replace(".hisat3n_dna.all_reads.name_sort.bam", "")
            bam_path = f"{sample_id}.hisat3n_dna.all_reads.name_sort.bam"
            output_prefix = f"{sample_id}.hisat3n_dna.all_reads"

            call_chromatin_contacts(
                bam_path=bam_path,
                contact_prefix=output_prefix,
                save_raw=False,
                save_hic_format=True
            )

        CODE

        #tar up the all_reads.contact_stats.csv files
        tar -zcvf ~{plate_id}.chromatin_contact_stats.tar.gz *.hisat3n_dna.all_reads.contact_stats.csv
        #tar up the .hisat3n_dna.all_reads.dedup_contacts.tsv files
        tar -zcvf ~{plate_id}.hisat3n_dna.all_reads.dedup_contacts.tar.gz *.hisat3n_dna.all_reads.dedup_contacts.tsv.gz
        #tar up the .hisat3n_dna.all_reads.3C.contact.tsv.gz files
        tar -zcvf ~{plate_id}.hisat3n_dna.all_reads.3C.contact.tar.gz *.hisat3n_dna.all_reads.3C.contact.tsv.gz
    >>>
    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
    }
    output {
        File chromatin_contact_stats = "~{plate_id}.chromatin_contact_stats.tar.gz"
        File all_reads_dedup_contacts = "~{plate_id}.hisat3n_dna.all_reads.dedup_contacts.tar.gz"
        File all_reads_3C_contacts = "~{plate_id}.hisat3n_dna.all_reads.3C.contact.tar.gz"
    }
}

task unique_reads_allc_and_cgn_extraction {
    input {
        File bam_and_index_tar
        File genome_fa
        String plate_id
        Int num_upstr_bases
        Int num_downstr_bases
        Int compress_level
        File chromosome_sizes

        Int disk_size = 200
        Int mem_size = 20
        String genome_base = basename(genome_fa)
        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
        Int preemptible_tries = 3
        Int cpu = 8
    }
    command <<<
        set -euo pipefail

        # unzip files
        tar -xf ~{bam_and_index_tar}
        rm ~{bam_and_index_tar}

        mkdir reference
        cp ~{genome_fa} reference
        cd reference

        # index the fasta
        echo "Indexing FASTA"
        samtools faidx *.fa
        cd ../output_bams

        echo "Starting allcools"
        bam_files=($(ls | grep "\.hisat3n_dna.all_reads.deduped.bam$"))
        echo ${bam_files[@]}
        for file in "${bam_files[@]}"; do
          sample_id=$(basename "$file" ".hisat3n_dna.all_reads.deduped.bam")
          /opt/conda/bin/allcools bam-to-allc \
          --bam_path "$file" \
          --reference_fasta /cromwell_root/reference/~{genome_base} \
          --output_path "${sample_id}.allc.tsv.gz" \
          --num_upstr_bases ~{num_upstr_bases} \
          --num_downstr_bases ~{num_downstr_bases} \
          --compress_level ~{compress_level} \
          --save_count_df \
          --convert_bam_strandness
        done
        echo "Zipping files"

        tar -zcvf ../~{plate_id}.allc.tsv.tar.gz *.allc.tsv.gz
        tar -zcvf ../~{plate_id}.allc.tbi.tar.gz *.allc.tsv.gz.tbi
        tar -zcvf ~{plate_id}.allc.count.tar.gz *.allc.tsv.gz.count.csv

        cd ../
        tar -xf ~{plate_id}.allc.tsv.tar.gz
        tar -xf ~{plate_id}.allc.tbi.tar.gz

        # prefix="allc-{mcg_context}/{cell_id}"
        if [ ~{num_upstr_bases} -eq 0 ]; then
           mcg_context=CGN
        else
           mcg_context=HCGN
        fi
        # create output dir
        mkdir /cromwell_root/allc-${mcg_context}
        outputdir=/cromwell_root/allc-${mcg_context}

        for gzfile in *.allc.tsv.gz
        do
             name=`echo $gzfile | cut -d. -f1`
             echo $name
             allcools extract-allc --strandness merge --allc_path $gzfile \
             --output_prefix $outputdir/$name \
             --mc_contexts ${mcg_context} \
             --chrom_size_path ~{chromosome_sizes}
        done

        mv output_bams/~{plate_id}.allc.count.tar.gz /cromwell_root

        cd /cromwell_root
        tar -zcvf ~{plate_id}.extract-allc.tar.gz $outputdir/*.gz
        tar -zcvf ~{plate_id}.extract-allc_tbi.tar.gz $outputdir/*.tbi

    >>>

    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
    }
    output {
        File allc = "~{plate_id}.allc.tsv.tar.gz"
        File tbi = "~{plate_id}.allc.tbi.tar.gz"
        File allc_uniq_reads_stats = "~{plate_id}.allc.count.tar.gz"
        File extract_allc_output_allc_tar = "~{plate_id}.extract-allc.tar.gz"
        File extract_allc_output_tbi_tar = "~{plate_id}.extract-allc_tbi.tar.gz"
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

        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
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

        python3 <<CODE
        from cemba_data.hisat3n import *
        snm3c_summary()
        CODE

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