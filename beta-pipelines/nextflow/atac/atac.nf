process GetNumSplits {
    publishDir "/mnt/disks/lk-disk/nextflow_test/"
    container = "ubuntu@sha256:2e863c44b718727c860746568e1d54afd13b2fa71b160f5cd9058fc436217b30"

    input:
    val nthreads
    val mem_size
    val cpu_platform

    output:
    path 'ranks_per_node.txt'
	
    script:
    """
    set -euo pipefail
    echo "Get number of splits for bwa-mem2"
    echo "#############################################"
    echo "Machine specs for bwa-mem2 task"
    echo "#############################################"

    num_nodes=1
    lscpu > compute_config

    num_cpus_per_node=\$(grep -E '^CPU\\(s\\)' compute_config | awk '{print \$2}')
    num_sockets=\$(grep -E '^Socket' compute_config | awk '{print \$2}')
    num_numa=\$(grep '^NUMA node(s)' compute_config | awk '{print \$3}')
    threads_per_core=\$(grep -E '^Thread' compute_config | awk '{print \$4}')

    num_cpus_all_node=\$((num_cpus_per_node * num_nodes))

    echo "Number of threads: " \$num_cpus_per_node
    echo "Number of sockets: " \$num_sockets
    echo "Number of NUMA domains: " \$num_numa
    echo "Number of threads per core: " \$threads_per_core
    echo "Number of CPUs: " \$num_cpus_all_node

    num_physical_cores_all_nodes=\$((num_cpus_all_node / threads_per_core))
    num_physical_cores_per_nodes=\$((num_cpus_per_node / threads_per_core))
    num_physical_cores_per_socket=\$((num_physical_cores_all_nodes / num_sockets))
    num_physical_cores_per_numa=\$((num_physical_cores_all_nodes / num_numa))

    echo "Number physical cores: " \$num_physical_cores_per_nodes
    echo "Number physical cores per socket: " \$num_physical_cores_per_socket
    echo "Number physical cores per numa: " \$num_physical_cores_per_numa

    th=\$((num_physical_cores_per_numa / 2))
    if [[ \$th -le 10 ]]; then
        th=\$num_physical_cores_per_numa
    fi

    while [[ \$num_physical_cores_per_nodes -gt \$th ]]; do
        num_physical_cores_per_nodes=\$((num_physical_cores_per_nodes / 2))
    done

    num_physical_cores_per_rank=\$num_physical_cores_per_nodes
    total_num_ranks=\$((num_physical_cores_all_nodes / num_physical_cores_per_rank))
    ranks_per_node=\$((total_num_ranks / num_nodes))

    echo "Number of MPI ranks: " \$total_num_ranks
    echo "Number of cores per MPI rank: " \$num_physical_cores_per_nodes
    echo "#############################################"
    echo "Please set bam_size such that fastqprocess creates \$total_num_ranks splits of input fastq files"
    echo "#############################################"

    echo \$total_num_ranks > total_num_ranks.txt
    echo \$ranks_per_node >> ranks_per_node.txt
    """
}

params.nthreads = 8
params.mem_size = "32GB"
params.cpu_platform = "x86_64"

workflow {
  result = GetNumSplits(params.nthreads, params.mem_size, params.cpu_platform)
  // Read and transform the output file into an integer
  ranks_per_node = result.ranks_per_node.text.trim().toInteger()
  println "Ranks per node: ${ranks_per_node}"

}