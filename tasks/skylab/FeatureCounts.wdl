version 1.0

task CountAlignments {

  input {
    File input_bam
    File annotation_gtf

    #runtime values
    String docker = "quay.io/humancellatlas/snss2-featurecount:0.1.0"
    Int machine_mem_mb = 8250
    Int cpu = 1
    Int disk = ceil(size(input_bam, "Gi") * 2) + 10
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

  command {
    set -e

   # counting the introns
   featureCounts -M -p ~{input_bam} \
      -a ~{annotation_gtf} \
      -F GTF \
      -o introns.counts \
      --minOverlap 3  \
      -t intron  \
      -g gene_id

   # create a new input bam where the alignemnts crossing intron-exon junctions are removed
   python /tools/remove-reads-on-junctions.py --input-gtf  ~{annotation_gtf} \
    --input-bam ~{input_bam}  --output-bam input.nojunc.bam

   # counting the exons
   featureCounts -M -p input.nojunc.bam \
    -a ~{annotation_gtf} -F GTF -o exons.counts  \
    --minOverlap 1 -t exon  -g gene_id 

  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File intron_counts_out = "introns.counts"
    File intron_counts_out_summary = "introns.counts.summary"
    File exon_counts_out = "exons.counts"
    File exon_counts_out_summary = "exons.counts.summary"
  }
}
