version 1.0

struct References {
  File genome_fa
  File annotation_gtf
}

task BuildStarSingleNucleus {
  input {
    # GTF annotation version refers to the version (GENCODE) or release (NCBI) listed in the GTF
    String gtf_annotation_version
    # Genome source can be NCBI or GENCODE
    String genome_source
    # Genome build is the assembly accession (NCBI) or version (GENCODE)
    String genome_build
    # Organism can be Macaque, Mouse, Human, etc.
    String organism
    File genome_fa
    File annotation_gtf
    File biotypes
    Int disk = 100
  }
  meta {
    description: "Modify GTF files and build reference index files for STAR aligner"
  }
  String ref_name = "star2.7.10a-~{organism}-~{genome_source}-build-~{genome_build}-~{gtf_annotation_version}"
  String star_index_name = "modified_~{ref_name}.tar"
  String annotation_gtf_modified = "modified_v~{gtf_annotation_version}.annotation.gtf"
  String annotation_gtf_introns = "introns_modified_v~{gtf_annotation_version}.annotation.gtf"

  command <<<
    # Check that input GTF files contain input genome source, genome build version, and annotation version
    if head -10 ~{annotation_gtf} | grep -qi ~{genome_build}
    then
        echo Genome version found in the GTF file
    else
        echo Error: Input genome version does not match version in GTF file
        exit 1;
    fi
    # Check that GTF file contains correct build source info in the first 10 lines of the GTF
    if head -10 ~{annotation_gtf} | grep -qi ~{genome_build}
    then
        echo Source of genome build identified in the GTF file
    else
        echo Error: Source of genome build not identified in the GTF file
        exit 1;
    fi

    set -eo pipefail

    python3 /script/modify_gtf.py  \
    --input-gtf ~{annotation_gtf} \
    --output-gtf ~{annotation_gtf_modified} \
    --biotypes ~{biotypes}

    python3  /script/add-introns-to-gtf.py   --input-gtf ~{annotation_gtf_modified}  --output-gtf ~{annotation_gtf_introns}

    mkdir star
    STAR --runMode genomeGenerate \
    --genomeDir star \
    --genomeFastaFiles ~{genome_fa} \
    --sjdbGTFfile ~{annotation_gtf_modified} \
    --sjdbOverhang 100 \
    --runThreadN 16

    tar -cvf ~{star_index_name} star

  >>>

  output {
    File star_index = star_index_name
    File annotation_gtf_modified_introns = annotation_gtf_introns
    File modified_annotation_gtf = annotation_gtf_modified
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/build-indices:1.0.0-2.7.10a-1671490724"
    memory: "50 GiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu:"16"
  }
}



workflow BuildIndices {
  input {
    # Genome source can be NCBI or GENCODE
    String genome_source
    # GTF annotation version refers to the version or release listed in the GTF
    String gtf_annotation_version
    # Genome build is the assembly accession (NCBI) or version (GENCODE)
    String genome_build
    File annotations_gtf
    File genome_fa
    String gtf_annotation_version
    File biotypes
  }

  # version of this pipeline
  String pipeline_version = "2.1.0"

  parameter_meta {
    annotations_gtf: "the annotation file"
    genome_fa: "the fasta file"
    biotypes: "gene_biotype attributes to include in the gtf file"
  }

  call BuildStarSingleNucleus {
    input:
      gtf_annotation_version = gtf_annotation_version,
      genome_fa = genome_fa,
      annotation_gtf = annotations_gtf,
      biotypes = biotypes,
      genome_build = genome_build,
      genome_source = genome_source
  }

  output {

    File snSS2_star_index = BuildStarSingleNucleus.star_index
    String pipeline_version_out = "BuildIndices_v~{pipeline_version}"
    File snSS2_annotation_gtf_introns = BuildStarSingleNucleus.annotation_gtf_modified_introns
    File snSS2_annotation_gtf_modified = BuildStarSingleNucleus.modified_annotation_gtf
  }
}
