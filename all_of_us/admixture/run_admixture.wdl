version 1.0

# For more information on the files and the contents, see: http://dalexander.github.io/admixture/admixture-manual.pdf
workflow run_admixture {
    input {
        File bed
        File bim
        File fam
    }

    String pipeline_version="aou_9.0.2"

    call run_admixture {
        input:
            bed=bed,
            bim=bim,
            fam=fam
    }

    output {
        File admixture_Q = run_admixture.Q
        File admixture_P = run_admixture.P
    }
}
task run_admixture {
    input {
        File bed
        File bim
        File fam
        Int? K_in
        Int? num_cpus_in
        Int mem_gb = 120
    }
    Int K = select_first([K_in, 6])
    Int num_cpus = select_first([num_cpus_in, 4])
    String basename = basename(bed, ".bed")

    command <<<
        set -e
        /app/bin/admixture ~{bed} ~{K} -j~{num_cpus}

        ls -la

    >>>

    output {
        File Q = "~{basename}.~{K}.Q"
        File P = "~{basename}.~{K}.P"
    }

    runtime {
        docker: "mussmann/admixpipe:3.0"
        memory: mem_gb + " GB" # Was 31 GB originally, increased for local ancestry
        cpu: "~{num_cpus}"
        disks: "local-disk 500 HDD"
    }
}