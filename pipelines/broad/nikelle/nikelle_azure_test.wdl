version 1.0

workflow Nikelle_azure_test {
    input {
   # File atac_whitelist = "gs://broad-gotc-test-storage/Multiome/input/737K-arc-v1_atac.txt"
   # File gex_whitelist = "gs://broad-gotc-test-storage/Multiome/input/737K-arc-v1_gex.txt"
    String cloud_provider
    File atac_whitelist
    File gex_whitelist
    }

String atac_whitelist_filename = "737K-arc-v1_atac.txt"
String gex_whitelist_filename = "737K-arc-v1_gex.txt"

String gcr_prefix = "gs://broad-gotc-test-storage/Multiome/input/"
String acr_prefix = "https://pdt656435155bfc8e.blob.core.windows.net/inputs/multiome/"

# choose docker prefix based on cloud provider
String cloud_prefix = if cloud_provider == "gcp" then gcr_prefix else acr_prefix

    call catThePaths as catThePaths_atac {
      input:
         whitelist = cloud_prefix + atac_whitelist_filename,
    }
    call catThePaths as catThePaths_gex {
      input:
          whitelist = cloud_prefix + gex_whitelist_filename,
    }

      output {
          File outputFile = catThePaths_atac.outPutFile
          File outputFile = catThePaths_gex.outPutFile
       }
}



task catThePaths {
    input{
        File whitelist
        String docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    }


    command <<<
        set -e
        set -o pipefail

        gsutil cat ~{whitelist} > whitelist.txt
    >>>

    runtime {
        docker: docker
        cpu : 1
        memory : "50 MiB"
        disks : "local-disk 10 HDD"
    }

    output {
        String outPutFile = read_string("whitelist.txt")
    }
}
