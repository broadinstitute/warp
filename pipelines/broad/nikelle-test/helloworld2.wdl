workflow HelloWorld2 {

    String s = "Hello World!"

    call echo {
        input: s = s
    }
}

task echo {

    String s

    command {
        picard.jar ${s}!
    }
    runtime {
        docker: "ubuntu"
        cpu: 1
        disks: "local-disk 1 HDD"
        preemptible: 1
        yo: "sup"
    }
    output {
        String out = read_string(stdout())
    }
}