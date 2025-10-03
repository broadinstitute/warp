version 1.0

task CountAlignments {

  input {
    Array[File] aligned_bam_inputs
    Array[String] input_ids
    File annotation_gtf

    #runtime values
    String docker = "us.gcr.io/broad-gotc-prod/subread:1.0.0-2.0.1-1689097353"
    Int machine_mem_mb = 8250
    Int cpu = 1
    Int disk = ceil(size(aligned_bam_inputs,"Gi")*2) + 10
    Int preemptible = 3
  }

  meta {
    description: "Counts the exonic and intronic reads from a bamfile using featureCounts."
  }

  parameter_meta {
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command <<<
    set -e

    declare -a bam_file_paths=(~{sep=' ' aligned_bam_inputs})
    declare -a output_prefix=(~{sep=' ' input_ids})

    # move the bam files to get around the issue of long file paths causing a segfault in featureCounts
    declare -a bam_files
    for filepath in ${bam_file_paths[@]};
      do
        filename="$(basename $filepath)"
        mv $filepath $filename
        bam_files+=("$filename")
      done;

    for (( i=0; i<${#bam_files[@]}; ++i));
      do
        # counting the introns
        featureCounts -M -p "${bam_files[$i]}" \
          -a ~{annotation_gtf} \
          -F GTF \
          -o "${output_prefix[$i]}.introns.counts" \
          --minOverlap 3  \
          -t intron  \
          -g gene_id

       # create a new input bam where the alignemnts crossing intron-exon junctions are removed
       python3 /usr/gitc/remove-reads-on-junctions.py --input-gtf  ~{annotation_gtf} \
        --input-bam "${bam_files[$i]}"  --output-bam "${output_prefix[$i]}.input.nojunc.bam"

       # counting the exons
       featureCounts -M -p "${output_prefix[$i]}.input.nojunc.bam" \
        -a ~{annotation_gtf} -F GTF -o "${output_prefix[$i]}.exons.counts"  \
        --minOverlap 1 -t exon  -g gene_id
      done;
  >>>

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu: cpu
    preemptible: preemptible
  }

  output {
    Array[File] intron_counts_out = glob("*.introns.counts")
    Array[File] intron_counts_out_summary = glob("*introns.counts.summary")
    Array[File] exon_counts_out = glob("*exons.counts")
    Array[File] exon_counts_out_summary = glob("*.exons.counts.summary")
  }
}
