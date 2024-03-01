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
        File r1_trimmed_fq_tar
        File r2_trimmed_fq_tar

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
    }

    # version of the pipeline
    String pipeline_version = "3.0.0"

    call Hisat_3n_pair_end_mapping_dna_mode {
        input:
            r1_trimmed_tar = r1_trimmed_fq_tar,
            r2_trimmed_tar = r2_trimmed_fq_tar,
            tarred_index_files = tarred_index_files,
            genome_fa = genome_fa,
            chromosome_sizes = chromosome_sizes,
            plate_id = plate_id
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

        String docker = "us.gcr.io/broad-gotc-prod/hisat3n:2.0.0-2.2.1-1708565445"
        Int disk_size = 1000
        Int mem_size = 64
        Int preemptible_tries = 3
        Int cpu = 48
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

        echo "starting hisat"

        task() {
            sample_id=$(basename "$file" "-R1_trimmed.fq.gz")
            hisat-3n /cromwell_root/$genome_fa_basename \
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
        }

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

        # Wait for all background jobs to finish before continuing
        wait

        echo "done hisat"

        echo "tarring up the outputs"
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