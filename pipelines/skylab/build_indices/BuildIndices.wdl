version 1.0

struct References {
  File genome_fa
  File annotation_gtf
}

task BuildStar {
  input {
    String gtf_version = "104"
    String organism = "Mmul10_SIV"
    References references
  }

  meta {
    description: "build reference index files for STAR aligner"
  }

  String ref_name = "star_primary_gencode_~{organism}_v~{gtf_version}"
  String star_index_name = "~{ref_name}.tar"

  command <<<
    set -eo pipefail

    mkdir star
    STAR --runMode genomeGenerate \
    --genomeDir star \
    --genomeFastaFiles ~{references.genome_fa} \
    --sjdbGTFfile ~{references.annotation_gtf} \
    --sjdbOverhang 100 \
    --runThreadN 16

    tar -cvf ~{star_index_name} star
  >>>

  output {
    File star_index = star_index_name
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-star:v2.7.9a"
    memory: "50 GiB"
    disks :"local-disk 100 HDD"
    cpu:"16"
  }
}

task BuildStarSingleNucleus {
  input {
    String gtf_version = "104"
    String organism = "Mmul10_SIV"
    String organism_prefix = "Mmul10_SIV"
    References references
    File biotypes
  }

  meta {
    description: "Modify gtf files and build reference index files for STAR aligner"
  }
  String ref_name = "star_primary_gencode_~{organism}_v~{gtf_version}"
  String star_index_name = "modified_~{ref_name}.tar"
  String annotation_gtf_modified = "modified_gencode.v~{gtf_version}.primary_assembly.annotation.gtf"
  String annotation_gtf_introns = "introns_modified_gencode.v~{gtf_version}.primary_assembly.annotation.gtf"

  command <<<
    set -eo pipefail

    python3 /script/modify_gtf.py  \
    --input-gtf ~{references.annotation_gtf} \
    --output-gtf ~{annotation_gtf_modified} \
    --biotypes ~{biotypes}

    python3  /script/add-introns-to-gtf.py   --input-gtf ~{annotation_gtf_modified}  --output-gtf ~{annotation_gtf_introns}

    mkdir star
    STAR --runMode genomeGenerate \
    --genomeDir star \
    --genomeFastaFiles ~{references.genome_fa} \
    --sjdbGTFfile ~{annotation_gtf_modified} \
    --sjdbOverhang 100 \
    --runThreadN 16

    tar -cvf ~{star_index_name} star

  >>>

  output {
    File star_index = star_index_name
    File annotation_gtf_modified_introns = annotation_gtf_introns
  }

  runtime {
    docker: "quay.io/humancellatlas/snss2-indices:1.2.0-scorch-rc1"
    memory: "50 GiB"
    disks :"local-disk 100 HDD"
    cpu:"16"
  }
}


workflow BuildIndices {
  input {
    File annotations_gtf
    File genome_fa
    String gtf_version = "104"
    String organism = "Mmul10_SIV"
    String organism_prefix = "Mmul10_SIV"
    File biotypes
  }

  # version of this pipeline

  String pipeline_version = "1.0.0"

  parameter_meta {
    gtf_version: "the actual number of gencode, ex.  27"
    organism: "Either 'human' or 'mouse'"
    organism_prefix: "Either 'h' or 'm'"
    biotypes: "gene_biotype attributes to include in the gtf file"
  }

  References references = object {
                            genome_fa: genome_fa,
                            annotation_gtf: annotations_gtf
                          }

  call BuildStar {
    input:
      gtf_version = gtf_version,
      organism = organism,
      references = references
  }

  call BuildStarSingleNucleus {
    input:
      gtf_version = gtf_version,
      organism = organism,
      organism_prefix = organism_prefix,
      references = references,
      biotypes = biotypes
  }

  output {
    File star_index = BuildStar.star_index
    File snSS2_star_index = BuildStarSingleNucleus.star_index
    String pipeline_version_out = "BuildIndices_v~{pipeline_version}"
    File snSS2_annotation_gtf_introns = BuildStarSingleNucleus.annotation_gtf_modified_introns
  }
}
