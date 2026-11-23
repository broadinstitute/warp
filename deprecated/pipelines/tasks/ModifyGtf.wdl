version 1.0

task ReplaceGeneNameWithGeneID {
  input {
    File original_gtf

    #runtime values
    String docker = "quay.io/humancellatlas/modify-gtf:0.1.0"
    Int machine_mem_mb = 8250
    Int cpu = 1
    Int disk = ceil(size(original_gtf, "Gi") * 2) + 10
    Int preemptible = 3
    String modified_gtf_location = "gene_id_as_gene_name.gtf.gz"
  }

  meta {
    description: "Modifies the gene_name field in a gtf to contain the values of gene_id instead."
  }

  parameter_meta {
    original_gtf: "The gtf to modify."
    docker: "(optional) the docker image containing the runtime environment for this task"
    modified_gtf_location: "(optional) the name to save the modified gtf file under"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -e

    if file --mime-type "${original_gtf}" | grep  gzip; then
      gunzip -c  "${original_gtf}" > input.gtf
    else
      mv "${original_gtf}"  input.gtf
    fi

    SetGeneNameToId.py \
      --in-gtf-file input.gtf \
      --out-gtf-file temp.gtf

    gzip -c temp.gtf > "${modified_gtf_location}"
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File modified_gtf = modified_gtf_location
  }
}
