version 1.0

workflow Methylome {
    meta {
        description: ""
        allowNestedInputs: true
    }

    input {
      File fastq_input_read1
      File fastq_input_read2
      File random_primer_indexes
    }


    call Demultiplex {
      input:
        fastq_input_read1 = fastq_input_read1,
        fastq_input_read2 = fastq_input_read2,
        random_primer_indexes = random_primer_indexes
    }

    output {
        Array[File] output_fastqs = Demultiplex.output_fastqs
    }
}

task Demultiplex {
    input {
        File fastq_input_read1
        File fastq_input_read2
        File random_primer_indexes

        String docker_image = "ekiernan/yap_hisat:v4"
        Int disk_size = 50
        Int mem_size = 10
}
    String base_name = basename(fastq_input_read1, ".fastq.gz")

command <<<
    set -euo pipefail

    /opt/conda/bin/cutadapt -Z -e 0.01 --no-indels -g file:~{random_primer_indexes} -o ~{base_name}.r1.{name}.output.fq.gz -p ~{base_name}.r2.{name}.output.fq.gz ~{fastq_input_read1} ~{fastq_input_read2} > stats.txt
>>>

runtime {
    docker: docker_image
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
    memory: "${mem_size} GiB"
}


output {
    Array[File] output_fastqs = glob("*.fq.gz")
}
}
