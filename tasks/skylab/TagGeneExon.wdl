version 1.0

task TagGeneExon {
  input {
    File annotations_gtf
    File bam_input

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-dropseqtools:v0.2.2-1.13"
    Int machine_mem_mb = 8250
    Int cpu = 1
    Int disk = ceil((size(bam_input, "Gi") + size(annotations_gtf, "Gi")) * 3) + 20
    Int preemptible = 3
  }

  meta {
    description: "Tags any read in bam_input that overlaps an intron or exon interval with the gene that those interals correspond to."
  }

  parameter_meta {
    annotations_gtf: "GTF annotation file for the species that bam input is derived from. Each record must have a gene_name and transcript_name in addition to a gene_id and transcript_id, no white space at the end of any record and must be in gtf format."
    bam_input: "Aligned bam file."
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

 command {
    set -e

    TagReadWithGeneExon \
      INPUT=${bam_input} \
      OUTPUT=bam_with_gene_exon.bam \
      SUMMARY=gene_exon_tag_summary.log \
      TAG=GE \
      ANNOTATIONS_FILE=${annotations_gtf}
  }

  # Larger genomes (mouse-human) require a 7.5gb instance; single-organism genomes work with 3.75gb
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File bam_output = "bam_with_gene_exon.bam"
    File log = "gene_exon_tag_summary.log"
  }
} 


task TagReadWithGeneFunction {
  input {
    File annotations_gtf
    File bam_input

    String gene_name_tag = "gn"
    String gene_strand_tag = "gs"
    String gene_function_tag = "gf"

    String use_strand_info = "true"

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-dropseqtools:2.3.0"
    Int machine_mem_mb = 8250
    Int cpu = 1
    Int disk = ceil((size(bam_input, "Gi") + size(annotations_gtf, "Gi")) * 3) + 20
    Int preemptible = 3
  }

  meta {
    description: "Tags any read in bam_input that overlaps an intron or exon interval with the gene that those interals correspond to."
  }

  parameter_meta {
    annotations_gtf: "GTF annotation file for the species that bam input is derived from. Each record must have a gene_name and transcript_name in addition to a gene_id and transcript_id, no white space at the end of any record and must be in gtf format."
    bam_input: "Aligned bam file."
    gene_name_tag: "the tag used to denote gene name in the bam (default: gn)"
    gene_strand_tag: "the tag used to denote gene strand in the bam (default: gs)"
    gene_function_tag: "the tag used to denote gene function (INTRONIC, EXONIC, ...) in the output bam (default: gf)"

    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

 command {
    set -e

    TagReadWithGeneFunction \
      INPUT=${bam_input} \
      OUTPUT=bam_with_gene_exon.bam \
      GENE_NAME_TAG=${gene_name_tag} \
      GENE_STRAND_TAG=${gene_strand_tag} \
      GENE_FUNCTION_TAG=${gene_function_tag} \
      SUMMARY=gene_exon_tag_summary.log \
      ANNOTATIONS_FILE=${annotations_gtf} \
      USE_STRAND_INFO=${use_strand_info}
  }

  # Larger genomes (mouse-human) require a 7.5gb instance; single-organism genomes work with 3.75gb
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File bam_output = "bam_with_gene_exon.bam"
    File log = "gene_exon_tag_summary.log"
  }
} 
