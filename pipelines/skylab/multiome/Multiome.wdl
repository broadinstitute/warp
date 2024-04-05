version 1.0

workflow Multiome {

    input {
        Boolean run_cellbender
    }

    # Call CellBender
    if (run_cellbender) {
        call print_cellbender {
            input:
                run_cellbender = run_cellbender
        }
    }

    output {
        File? run_cellbender_out1 = print_cellbender.out1
        File? run_cellbender_out2 = print_cellbender.out2
    }
}

task print_cellbender {
    input{
        Boolean run_cellbender
    }
    command <<<
        echo "~{run_cellbender}" > 1.txt
        echo ~{run_cellbender} > 2.txt
    >>>
    output {
        File out1 = "1.txt"
        File out2 = "2.txt"
    }
    runtime {
        docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
        disks: "local-disk 10 HDD"
        memory: "4 GB"
    }
}


