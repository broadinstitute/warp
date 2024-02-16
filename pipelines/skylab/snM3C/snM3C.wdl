version 1.0

workflow snM3C {

    input {
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

        File r1_trimmed_tar
        File r2_trimmed_tar
    }

    # version of the pipeline
    String pipeline_version = "2.0.0"

    call Hisat_3n_pair_end_mapping_dna_mode {
        input:
            r1_trimmed_tar = r1_trimmed_tar,
            r2_trimmed_tar = r2_trimmed_tar,
            tarred_index_files = tarred_index_files,
            genome_fa = genome_fa,
            chromosome_sizes = chromosome_sizes,
            plate_id = plate_id
    }

    output {
        File hisat3n_paired_end_bam_tar = Hisat_3n_pair_end_mapping_dna_mode.hisat3n_paired_end_bam_tar
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
        Int disk_size = 2000
        Int mem_size = 512
        Int preemptible_tries = 3
        Int cpu = 128
        String cpuPlatform = "Intel Ice Lake"
    }
    command <<<
        set -euo pipefail

        # check genomic reference version and print to output txt file
        STRING=~{genome_fa}
        BASE=$(basename $STRING .fa)

        echo "The reference is $BASE" > ~{plate_id}.reference_version.txt

        # untar the index files
        echo "Untarring the index files"
        date
        tar -zxvf ~{tarred_index_files}
        rm ~{tarred_index_files}
        echo "tarring has finished"
        date

        echo "copying genome"
        date
        cp ~{genome_fa} .
        echo "done copying genome"
        date

        #get the basename of the genome_fa file
        genome_fa_basename=$(basename ~{genome_fa} .fa)
        echo "samtools faidx $genome_fa_basename.fa"
        samtools faidx $genome_fa_basename.fa

        # untar the demultiplexed fastq files
        echo "Untarring the fastq files"
        date
        tar -zxvf ~{r1_trimmed_tar}
        tar -zxvf ~{r2_trimmed_tar}
        rm ~{r1_trimmed_tar}
        rm ~{r2_trimmed_tar}
        echo "done tarring the fastq files"
        date

        # define lists of r1 and r2 fq files
        R1_files=($(ls | grep "\-R1_trimmed.fq.gz"))
        R2_files=($(ls | grep "\-R2_trimmed.fq.gz"))

        # check to make sure these arrays are the same length
        #if [ ${#R1_files[@]} -ne ${#R2_files[@]} ]; then
        #  echo "The number of R1 and R2 files are not the same"
        #  exit 1
        #fi

        #turn the arrays into comma separated strings
        #R1_files_string=$(IFS=,; echo "${R1_files[*]}")
        #R2_files_string=$(IFS=,; echo "${R2_files[*]}")

        echo "starting hisat"
        date

        # Print the list of sample IDs
        #echo "List of Sample IDs:"
        #for id in "${sample_ids[@]}"; do
            #echo "$id"
        #done


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
            --threads 126 | samtools view -@ 10 -b -q 0 -o "${sample_id}.hisat3n_dna.unsort.bam" &
        done

        # Wait for all background jobs to finish before continuing
        wait

        echo "done hisat"
        date

        echo "tarring up the outputs"
        date
        # tar up the bam files and stats files
        tar -zcvf ~{plate_id}.hisat3n_paired_end_bam_files.tar.gz *.bam
        tar -zcvf ~{plate_id}.hisat3n_paired_end_stats_files.tar.gz *.hisat3n_dna_summary.txt
        echo "tarring up the outputs"
        date

    >>>
    runtime {
        docker: docker
        disks: "local-disk ${disk_size} HDD"
        cpu: cpu
        memory: "${mem_size} GiB"
        preemptible: preemptible_tries
        cpuPlatform: cpuPlatform
    }
    output {
        File hisat3n_paired_end_bam_tar = "~{plate_id}.hisat3n_paired_end_bam_files.tar.gz"
        File hisat3n_paired_end_stats_tar = "~{plate_id}.hisat3n_paired_end_stats_files.tar.gz"
        File reference_version = "~{plate_id}.reference_version.txt"
    }
}
