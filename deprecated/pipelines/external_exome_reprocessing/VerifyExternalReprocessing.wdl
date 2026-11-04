version 1.0

workflow VerifyExternalReprocessing {

    input {
        Boolean done
        String destination_cloud_path
    }

    call AssertTrue {
        input:
            done                   = done,
            destination_cloud_path = destination_cloud_path
    }

}

task AssertTrue {
    input {
        Boolean done
        String destination_cloud_path
    }

    command {
        echo "Verifying that ExternalReprocessing copy step completed..."
        if [! ~{done} ]; then  
            echo "ExternalReprocessing failed to copy to ${destination_cloud_path}" &>2
            exit 1
        fi

        exit 0
    }

    runtime {
        docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
        disks: "local-disk 10 HDD"
        memory: "3.5 GiB"
        preemptible: 3
    }
}