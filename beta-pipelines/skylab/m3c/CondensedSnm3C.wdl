version 1.0

workflow WDLized_snm3C {

    input {
        Array[File] fastq_input_read1
        Array[File] fastq_input_read2
        File random_primer_indexes
        String plate_id
        # mapping inputs
        File tarred_index_files
        File genome_fa
        File chromosome_sizes
    }

    call Demultiplexing {
        input:
            fastq_input_read1 = fastq_input_read1,
            fastq_input_read2 = fastq_input_read2,
            random_primer_indexes = random_primer_indexes,
            plate_id = plate_id
    }

    call Sort_and_trim_r1_and_r2 {
        input:
            tarred_demultiplexed_fastqs = Demultiplexing.tarred_demultiplexed_fastqs
    }

    call Hisat_3n_pair_end_mapping_dna_mode {
        input:
            r1_trimmed_tar = Sort_and_trim_r1_and_r2.r1_trimmed_fq_tar,
            r2_trimmed_tar = Sort_and_trim_r1_and_r2.r2_trimmed_fq_tar,
            tarred_index_files = tarred_index_files,
            genome_fa = genome_fa,
            chromosome_sizes = chromosome_sizes
    }

    call Separate_unmapped_reads {
        input:
            hisat3n_bam_tar = Hisat_3n_pair_end_mapping_dna_mode.hisat3n_paired_end_bam_tar
    }

    call Split_unmapped_reads {
        input:
            unmapped_fastq_tar = Separate_unmapped_reads.unmapped_fastq_tar
    }

    call Hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name {
        input:
            split_fq_tar = Split_unmapped_reads.split_fq_tar,
            tarred_index_files = tarred_index_files,
            genome_fa = genome_fa
    }

   call remove_overlap_read_parts {
       input:
           bam = Hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name.merge_sorted_bam_tar
   }

   # call merge_original_and_split_bam_and_sort_all_reads_by_name_and_position {
   #     input:
   #         unique_bam = separate_unmapped_reads.unique_bam,
   #         split_bam = remove_overlap_read_parts.remove_overlap_bam
   # }
#
   # call call_chromatin_contacts {
   #     input:
   #         bam = merge_original_and_split_bam_and_sort_all_reads_by_name_and_position.name_sorted_bam
   # }
#
   # call dedup_unique_bam_and_index_unique_bam {
   #     input:
   #         bam = merge_original_and_split_bam_and_sort_all_reads_by_name_and_position.position_sorted_bam
   # }
#
   # call unique_reads_allc {
   #     input:
   #         bam = dedup_unique_bam_and_index_unique_bam.dedup_bam,
   #         bai = dedup_unique_bam_and_index_unique_bam.dedup_bai
   # }
#
   # call unique_reads_cgn_extraction {
   #     input:
   #         allc = unique_reads_allc.allc,
   #         tbi = unique_reads_allc.tbi
   # }
#
   # call summary {
   #     input:
   #         trimmed_stats = Sort_and_trim_r1_and_r2.trim_stats,
   #         hisat3n_stats = Hisat_3n_pair_end_mapping_dna_mode.hisat3n_stats,
   #         r1_hisat3n_stats = hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name.r1_hisat3n_stats,
   #         r2_hisat3n_stats = hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name.r2_hisat3n_stats,
   #         dedup_stats = dedup_unique_bam_and_index_unique_bam.dedup_stats,
   #         chromatin_contact_stats = call_chromatin_contacts.chromatin_contact_stats,
   #         allc_uniq_reads_stats = unique_reads_allc.allc_uniq_reads_stats,
   #         unique_reads_cgn_extraction_tbi = unique_reads_cgn_extraction.unique_reads_cgn_extraction_tbi
   # }

    output {
        #File MappingSummary = summary.MappingSummary
        #File allcFiles = unique_reads_allc.
        #File allc_CGNFiles = unique_reads_cgn_extraction.
        #File bamFiles = merge_original_and_split_bam_and_sort_all_reads_by_name_and_position.name_sorted_bam
        #File UniqueAlign_cell_parser_picard_dedup = dedup_unique_bam_and_index_unique_bam.dedup_stats
        #File SplitReads_cell_parser_hisat_summary = "?"
        #File hicFiles = call_chromatin_contacts.chromatin_contact_stats
        File trimmed_stats = Sort_and_trim_r1_and_r2.trim_stats_tar
        File r1_trimmed_fq = Sort_and_trim_r1_and_r2.r1_trimmed_fq_tar
        File r2_trimmed_fq = Sort_and_trim_r1_and_r2.r2_trimmed_fq_tar
        File hisat3n_stats_tar = Hisat_3n_pair_end_mapping_dna_mode.hisat3n_paired_end_stats_tar
        File hisat3n_bam_tar = Hisat_3n_pair_end_mapping_dna_mode.hisat3n_paired_end_bam_tar
        File unique_bam_tar = Separate_unmapped_reads.unique_bam_tar
        File multi_bam_tar = Separate_unmapped_reads.multi_bam_tar
        File unmapped_fastq_tar = Separate_unmapped_reads.unmapped_fastq_tar
        File split_fq_tar = Split_unmapped_reads.split_fq_tar
        File merge_sorted_bam_tar = Hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name.merge_sorted_bam_tar
        File remove_overlap_read_parts_bam_tar = remove_overlap_read_parts.output_bam_tar
    }
}

task Demultiplexing {
  input {
    Array[File] fastq_input_read1
    Array[File] fastq_input_read2
    File random_primer_indexes
    String plate_id

    String docker_image = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
    Int disk_size = 50
    Int mem_size = 10
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

    # zip up all the output fq.gz files
    tar -zcvf ~{plate_id}.cutadapt_output_files.tar.gz *.fq.gz
  >>>

  runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
    memory: "${mem_size} GiB"
  }

  output {
    File tarred_demultiplexed_fastqs = "~{plate_id}.cutadapt_output_files.tar.gz"
    File stats = "~{plate_id}.stats.txt"
    }
}

task Sort_and_trim_r1_and_r2 {
    input {
        File tarred_demultiplexed_fastqs
        Int disk_size = 50
        Int mem_size = 10

        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"

    }
    command <<<
    set -euo pipefail

    # untar the demultiplexed fastqs
    tar -xf ~{tarred_demultiplexed_fastqs}

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
        -a R1Adapter='AGATCGGAAGAGCACACGTCTGAAC' \
        -A R2Adapter='AGATCGGAAGAGCGTCGTGTAGGGA' \
        --report=minimal \
        -O 6 \
        -q 20 \
        -u 10 \
        -u -10 \
        -U 10 \
        -U -10 \
        -Z \
        -m 30:30 \
        --pair-filter 'both' \
        -o ${sample_id}-R1_trimmed.fq.gz \
        -p ${sample_id}-R2_trimmed.fq.gz \
        ${sample_id}-R1_sorted.fq ${sample_id}-R2_sorted.fq \
        > ${sample_id}.trimmed.stats.txt
    done

    echo "Tarring up the trimmed files and stats files"

    tar -zcvf R1_trimmed_files.tar.gz *-R1_trimmed.fq.gz
    tar -zcvf R2_trimmed_files.tar.gz *-R2_trimmed.fq.gz
    tar -zcvf trimmed_stats_files.tar.gz *.trimmed.stats.txt
>>>
    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File r1_trimmed_fq_tar = "R1_trimmed_files.tar.gz"
        File r2_trimmed_fq_tar = "R2_trimmed_files.tar.gz"
        File trim_stats_tar = "trimmed_stats_files.tar.gz"
    }
}

task Hisat_3n_pair_end_mapping_dna_mode{
    input {
        File r1_trimmed_tar
        File r2_trimmed_tar
        File tarred_index_files
        File genome_fa
        File chromosome_sizes

        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
        Int disk_size = 100
        Int mem_size = 100
    }
    command <<<
        set -euo pipefail

        mkdir reference/
        mkdir fastq/

        cp ~{tarred_index_files} reference/
        cp ~{genome_fa} reference/
        cp ~{chromosome_sizes} reference/
        cp ~{r1_trimmed_tar} fastq/
        cp ~{r2_trimmed_tar} fastq/

        # untar the index files
        cd reference/
        echo "Untarring the index files"
        tar -zxvf ~{tarred_index_files}
        rm ~{tarred_index_files}
        samtools faidx hg38.fa

        # untar the demultiplexed fastq files
        cd ../fastq/
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
          hisat-3n /cromwell_root/reference/hg38 \
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
        tar -zcvf hisat3n_paired_end_bam_files.tar.gz *.bam
        tar -zcvf hisat3n_paired_end_stats_files.tar.gz *.hisat3n_dna_summary.txt

        mv hisat3n_paired_end_bam_files.tar.gz ../
        mv hisat3n_paired_end_stats_files.tar.gz ../

    >>>
    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File hisat3n_paired_end_bam_tar = "hisat3n_paired_end_bam_files.tar.gz"
        File hisat3n_paired_end_stats_tar = "hisat3n_paired_end_stats_files.tar.gz"
    }
}

task Separate_unmapped_reads {
    input {
        File hisat3n_bam_tar

        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
        Int disk_size = 50
        Int mem_size = 10

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
                qlen_cutoff=30
            )

        CODE

        # tar up the uniqe bams
        tar -zcvf hisat3n_paired_end_unique_bam_files.tar.gz *.hisat3n_dna.unique_aligned.bam

        # tar up the multi bams
        tar -zcvf hisat3n_paired_end_multi_bam_files.tar.gz *.hisat3n_dna.multi_aligned.bam

        # tar up the unmapped fastq files
        tar -zcvf hisat3n_paired_end_unmapped_fastq_files.tar.gz *.hisat3n_dna.unmapped.fastq

    >>>
    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File unique_bam_tar = "hisat3n_paired_end_unique_bam_files.tar.gz"
        File multi_bam_tar = "hisat3n_paired_end_multi_bam_files.tar.gz"
        File unmapped_fastq_tar = "hisat3n_paired_end_unmapped_fastq_files.tar.gz"
    }
}

task Split_unmapped_reads {
    input {
        File unmapped_fastq_tar
        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
        Int disk_size = 50
        Int mem_size = 10
    }
    command <<<

        set -euo pipefail

        # untar the unmapped fastq files
        tar -xf ~{unmapped_fastq_tar}
        rm ~{unmapped_fastq_tar}

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
            min_length=30
          )

        CODE

        # wait 15 seconds for the files to be written
        sleep 15

        # tar up the split fastq files
        tar -zcvf hisat3n_paired_end_split_fastq_files.tar.gz *.split_reads*.fastq


    >>>
    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File split_fq_tar = "hisat3n_paired_end_split_fastq_files.tar.gz"
    }
}

task Hisat_single_end_r1_r2_mapping_dna_mode_and_merge_sort_split_reads_by_name {
    input {
        File split_fq_tar
        File genome_fa
        File tarred_index_files

        String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
        Int disk_size = 80
        Int mem_size = 20
    }
    command <<<
        set -euo pipefail

        mkdir reference/

        cp ~{tarred_index_files} reference/
        cp ~{genome_fa} reference/

        # untar the tarred index files
        cd reference/
        tar -xvf ~{tarred_index_files}
        rm ~{tarred_index_files}
        samtools faidx hg38.fa

        # untar the unmapped fastq files
        tar -xvf ~{split_fq_tar}
        rm ~{split_fq_tar}

        # define lists of r1 and r2 fq files
        R1_files=($(ls | grep "\.hisat3n_dna.split_reads.R1.fastq"))
        R2_files=($(ls | grep "\.hisat3n_dna.split_reads.R2.fastq"))

        for file in "${R1_files[@]}"; do
          sample_id=$(basename "$file" ".hisat3n_dna.split_reads.R1.fastq")
          hisat-3n /cromwell_root/reference/hg38 \
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
         hisat-3n /cromwell_root/reference/hg38 \
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
       tar -zcvf ../hisat3n_dna.split_reads.name_sort.bam.tar.gz *.hisat3n_dna.split_reads.name_sort.bam

    >>>
    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }
    output {
        File merge_sorted_bam_tar = "hisat3n_dna.split_reads.name_sort.bam.tar.gz"
    }
}

task remove_overlap_read_parts {
   input {
       File bam
       String docker = "us.gcr.io/broad-gotc-prod/m3c-yap-hisat:1.0.0-2.2.1"
       Int disk_size = 80
       Int mem_size = 20  
   }

   command <<<
       set -euo pipefail
       # unzip bam file
       tar -xf ~{bam}
       rm ~{bam}

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
       tar -zcvf ../remove_overlap_read_parts.tar.gz *bam
       
   >>>
   runtime {
       docker: docker
       disks: "local-disk ${disk_size} HDD"
       cpu: 1
       memory: "${mem_size} GiB"
   }
   output {
        File output_bam_tar = "remove_overlap_read_parts.tar.gz"
    }

}

#task merge_original_and_split_bam_and_sort_all_reads_by_name_and_position {
#    input {
#        File unique_bam
#        File split_bam
#    }
#    command <<<
#    >>>
#    runtime {
#        docker: "fill_in"
#        disks: "local-disk ${disk_size} HDD"
#        cpu: 1
#        memory: "${mem_size} GiB"
#    }
#    output {
#        File name_sorted_bam = ""
#        File position_sorted_bam = ""
#    }
#}

#task call_chromatin_contacts {
#    input {
#        File bam
#    }
#    command <<<
#    >>>
#    runtime {
#        docker: "fill_in"
#        disks: "local-disk ${disk_size} HDD"
#        cpu: 1
#        memory: "${mem_size} GiB"
#    }
#    output {
#        File chromatin_contact_stats = ""
#    }
#}

#task dedup_unique_bam_and_index_unique_bam {
#    input {
#        File bam
#    }
#    command <<<
#    >>>
#    runtime {
#        docker: "fill_in"
#        disks: "local-disk ${disk_size} HDD"
#        cpu: 1
#        memory: "${mem_size} GiB"
#    }
#    output {
#        File dedup_bam = ""
#        File dedup_stats = ""
#        File dedup_bai = ""
#    }
#}

#task unique_reads_allc {
#    input {
#        File bam
#        File bai
#    }
#    command <<<
#    >>>
#    runtime {
#        docker: "fill_in"
#        disks: "local-disk ${disk_size} HDD"
#        cpu: 1
#        memory: "${mem_size} GiB"
#    }
#    output {
#        File allc = ""
#        File tbi = ""
#        File allc_uniq_reads_stats = ""
#    }
#}

#task unique_reads_cgn_extraction {
#    input {
#        File allc
#        File tbi
#    }
#    command <<<
#    >>>
#    runtime {
#        docker: "fill_in"
#        disks: "local-disk ${disk_size} HDD"
#        cpu: 1
#        memory: "${mem_size} GiB"
#    }
#    output {
#        File unique_reads_cgn_extraction_allc = ""
#        File unique_reads_cgn_extraction_tbi = ""
#    }
#}

#task summary {
#    input {
#        File trimmed_stats
#        File hisat3n_stats
#        File r1_hisat3n_stats
#        File r2_hisat3n_stats
#        File dedup_stats
#        File chromatin_contact_stats
#        File allc_uniq_reads_stats
#        File unique_reads_cgn_extraction_tbi
#    }
#    command <<<
#    >>>
#    runtime {
#        docker: "fill_in"
#        disks: "local-disk ${disk_size} HDD"
#        cpu: 1
#        memory: "${mem_size} GiB"
#    }
#    output {
#        File mapping_summary = "MappingSummary.csv.gz"
#    }
#}
