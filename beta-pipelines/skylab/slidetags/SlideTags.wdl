version 1.0

import "scripts/spatial-count.wdl" as SpatialCount
import "scripts/positioning.wdl" as Positioning

workflow SlideTags {

    String pipeline_version = "1.0.0"

    input {
        String id
        Array[String] fastq_paths
        Array[String] pucks
        Array[String] rna_paths
        String sb_path
        String docker = "us.gcr.io/broad-gotc-prod/slide-tags:1.1.0"
     }
    
    parameter_meta {
        fastq_paths: "Array of paths to spatial fastq files"
        pucks: "Array of paths to puck files"
        docker: "Docker image to use"
    }

    call SpatialCount.count as spatial_count {
        input:
            fastq_paths = fastq_paths,
            pucks = pucks,
            docker = docker
     }

    call Positioning.generate_positioning as positioning {
        input:
            rna_paths = rna_paths,
            sb_path = spatial_count.sb_counts,
            docker = docker
     }
    
}

