version 1.0

workflow Methylome {

    input {

        # mapping inputs
        File tarred_index_files
        File mapping_yaml
        File snakefile
        File chromosome_sizes
        File genome_fa
        File tarred_demultiplexed_fastqs = "gs://broad-gotc-test-storage/AD3C_BA17_2027_P1-1-B11-A1.cutadapt_output_files.tar.gz"
    }


    call Mapping {
        input:
            tarred_demultiplexed_fastqs = tarred_demultiplexed_fastqs,
            tarred_index_files = tarred_index_files,
            mapping_yaml = mapping_yaml,
            snakefile = snakefile,
            chromosome_sizes = chromosome_sizes,
            genome_fa = genome_fa
    }

    output {
        #File MappingSummary = Mapping.MappingSummary
        #Array[file] allcFiles = Mapping.allcFiles
        #Array[file] allc_CGNFiles = Mapping.allc_CGNFiles
        #Array[file] bamFiles = Mapping.bamFiles
        #Array[file] detail_statsFiles = Mapping.detail_statsFiles
        #Array[file] hicFiles = Mapping.hicFiles

    }
}



task Mapping {
    input {
        File tarred_demultiplexed_fastqs
        File tarred_index_files
        File mapping_yaml
        File snakefile
        File chromosome_sizes
        File genome_fa


        String docker_image = "nikellepetrillo/yap-hisat:v7"
        Int disk_size = 200
        Int mem_size = 500
    }

    command <<<
        set -euo pipefail

        #mkdir group0/
        #cd group0/
        #echo "pwd is"
        #pwd
        echo "call cutadapt"
        cutadapt --version

        mkdir group0/
        mkdir group0/fastq/
        mkdir group0/reference/

        cp ~{tarred_index_files} group0/reference/
        cp ~{chromosome_sizes} group0/reference/
        cp ~{genome_fa} group0/reference/
        cp ~{tarred_demultiplexed_fastqs} group0/fastq/
        cp ~{mapping_yaml} group0/
        cp ~{snakefile} group0/

        # untar the index files
        cd group0/reference/
        echo "Untarring the index files"
        tar -zxvf ~{tarred_index_files}
        rm ~{tarred_index_files}
        samtools faidx hg38.fa
        echo "The current working directory is (for the reference dir):"
        pwd
        echo "here is the ls command (for the reference dir):"
        ls
        echo "echo the path"
        echo $PATH



        # untar the demultiplexed fastq files
        cd ../fastq/
        echo "Untarring the demultiplexed fastq files"
        tar -zxvf ~{sep=' ' tarred_demultiplexed_fastqs}


        # run the snakemake command
        cd ../
        echo "The current working directory is  (the snakemake command is being run here:"
        pwd
        echo "here is the ls command (for the snakemake command):"
        ls

        /opt/conda/bin/snakemake --verbose --configfile mapping.yaml -j

        ls -l


    >>>

    runtime {
        docker: docker_image
        disks: "local-disk ${disk_size} HDD"
        cpu: 1
        memory: "${mem_size} GiB"
    }

    output {
        File MappingSummary = "/cromwell_root/group0/MappingSummary.csv.gz"
        Array[File] allcFiles = glob("/cromwell_root/group0/allc/*")
        Array[File] allc_CGNFiles = glob("/cromwell_root/group0/allc-CGN/*")
        Array[File] bamFiles = glob("/cromwell_root/group0/bam/*")
        Array[File] detail_statsFiles = glob("/cromwell_root/group0/detail_stats/*")
        Array[File] hicFiles = glob("/cromwell_root/group0/hic/*")






    }
}
