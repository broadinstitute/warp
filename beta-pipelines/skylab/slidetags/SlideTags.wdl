version 1.0

import "spatial-count/spatial-count.wdl" as SpatialCount

workflow SlideTags {

    String pipeline_version = "1.0.0"

    input {
        String id
        Array[String] fastq_paths
        Array[String] pucks
        Int mem_GiB = 64
        Int disk_GiB = 128
        String docker = "us-central1-docker.pkg.dev/velina-208320/docker-count/img:latest"
     }
    
    parameter_meta {
        fastq_paths: "Array of paths to spatial fastq files"
        pucks: "Array of paths to puck files"
        mem_GiB: "Memory in GiB to allocate to the task"
        disk_GiB: "Disk in GiB to allocate to the task"
        docker: "Docker image to use"
    }

    call SpatialCount.count as spatial_count {
        input:
            fastq_paths = fastq_paths,
            pucks = pucks,
            mem_GiB = mem_GiB,
            disk_GiB = disk_GiB,
            docker = docker
     }
    
}

