version 1.0

workflow LocalizeAndConcatStargazer {
    input {
        File list_of_gene_names
        File list_of_gene_fofns
    }
    
    String pipeline_version = "aou_9.0.0"
    Array[String] gene_names = read_lines(list_of_gene_names)
    Array[File] gene_fofns = read_lines(list_of_gene_fofns)
    Array[Pair[String, File]] task_pairs = zip(gene_names, gene_fofns)

    scatter(pair in task_pairs) {
        call LocalizeAndConcat {
            input:
                output_root = pair.left,
                fofn = pair.right
        }
    }

    output {
        Array[File] per_gene_results = LocalizeAndConcat.all_samples_output
    }
}

task LocalizeAndConcat {
    input {
        String output_root
        File fofn
    }

    command <<<
        set -euxo pipefail
        rm -rf per_sample
        mkdir per_sample
        cat ~{fofn} | gsutil -m cp -I per_sample/

        ls -1 per_sample > per_sample.list
        cd per_sample

        head -n1 $(head -n1 ../per_sample.list) > all_samples_~{output_root}.tsv

        while read output;
        do tail -n1 $output >> all_samples_~{output_root}.tsv;
        done < ../per_sample.list
    >>>

    runtime {
        cpu: 4
        memory: "6 GiB"
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
        disks: "local-disk 10 HDD"
    }

    output {
        File all_samples_output = "per_sample/all_samples_~{output_root}.tsv"
    }
}