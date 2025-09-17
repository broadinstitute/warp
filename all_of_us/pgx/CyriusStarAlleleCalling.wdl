version 1.0

# assume hg38

workflow CyriusStarAlleleCalling {
    input {
        File input_cram
        File input_cram_index
        String sample_name
        File interval_list
        File ref_dict
        File ref_fasta
        File ref_fasta_index
    }
    String pipeline_version = "aou_9.0.0"

    call PrintReads {
        input:
            input_bam = input_cram,
            input_bam_index = input_cram_index,
            sample_name = sample_name,
            interval_list = interval_list,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
    }

    call RunCyrius {
        input:
            input_cram = PrintReads.subset_bam,
            input_cram_index = PrintReads.subset_bam_index,
            sample_name = sample_name,
            ref_fasta = ref_fasta
    }

    output {
        File subset_cram = PrintReads.subset_bam
        File subset_cram_index = PrintReads.subset_bam_index
        File cyrius_output = RunCyrius.cyrius_output
        File cyrius_details = RunCyrius.cyrius_details
    }
}

task PrintReads {
    input {
        File input_bam
        File input_bam_index
        String sample_name
        File interval_list
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        Int preemptible_tries = 2
        String requester_pays_project
    }

    parameter_meta {
        input_bam: {
                       localization_optional: true
                   }
    }

    Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
    Int disk_size = ceil(ref_size) + 20

    command {
        /usr/gitc/gatk4/gatk --java-options "-Xms2000m -Xmx9000m" \
        PrintReads \
        -I ~{input_bam} \
        -R ~{ref_fasta} \
        --interval-padding 500 \
        -L ~{interval_list} \
        -O ~{sample_name}.cyp2d6_reads.bam 
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/gatk:1.3.0-4.2.6.1-1649964384"
        preemptible: preemptible_tries
        memory: "10000 MiB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
    }
    output {
        File subset_bam = "~{sample_name}.cyp2d6_reads.bam"
        File subset_bam_index = "~{sample_name}.cyp2d6_reads.bai"
    }
}

task RunCyrius {
    input {
        File input_cram
        File input_cram_index
        String sample_name
        File ref_fasta
    }

    command {
        echo ~{input_cram} > manifest.txt

        python /Cyrius/star_caller.py --manifest manifest.txt \
        --genome 38 \
        --prefix ~{sample_name} \
        --outDir ~{sample_name} \
        --threads 2 --reference ~{ref_fasta}
    }

    output {
        File cyrius_output = "~{sample_name}/~{sample_name}.tsv"
        File cyrius_details = "~{sample_name}/~{sample_name}.json"
    }

    runtime {
        memory: "3750 MiB"
        disks: "local-disk 100 HDD"
        bootDiskSizeGb: 15
        preemptible: 3
        docker: "us.gcr.io/broad-dsde-methods/cyrius@sha256:346586f32d644d26014163b7c8799b53b8ce501af558405e01f6b8456de0ca45"
    }
}