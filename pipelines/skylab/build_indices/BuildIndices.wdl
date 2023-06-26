version 1.0

workflow BuildIndices {
  input {
    # Genome source can be NCBI or GENCODE
    String genome_source
    # GTF annotation version refers to the version or release listed in the GTF
    String gtf_annotation_version
    # Genome build is the assembly accession (NCBI) or version (GENCODE)
    String genome_build
    String organism

    File annotations_gtf
    File genome_fa
    File biotypes
  }

  # version of this pipeline
  String pipeline_version = "2.1.2"


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
        genome_source = genome_source,
        organism = organism
    }
    call CalculateChromosomeSizes {
      input:
        genome_fa = genome_fa
    }
    call BuildBWAreference {
      input:
        genome_fa = genome_fa,
        chrom_sizes_file = CalculateChromosomeSizes.chrom_sizes,
        genome_source = genome_source,
        genome_build = genome_build,
        gtf_annotation_version = gtf_annotation_version,
        organism = organism
    }

  output {
    File snSS2_star_index = BuildStarSingleNucleus.star_index
    String pipeline_version_out = "BuildIndices_v~{pipeline_version}"
    File snSS2_annotation_gtf_introns = BuildStarSingleNucleus.annotation_gtf_modified_introns
    File snSS2_annotation_gtf_modified = BuildStarSingleNucleus.modified_annotation_gtf
    File reference_bundle = BuildBWAreference.reference_bundle
    File chromosome_sizes = CalculateChromosomeSizes.chrom_sizes
  }
}

task CalculateChromosomeSizes {
  input {
    File genome_fa
  }
  command <<<
    samtools faidx ~{genome_fa} | cut -f1,2 > chrom.sizes
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
    preemptible: 3
    memory: "3 GiB"
    cpu: "1"
    disks: "local-disk 50 HDD"
  }
  output {
    File chrom_sizes = "chrom.sizes"
  }
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
    if head -10 ~{annotation_gtf} | grep -qi ~{genome_source}
    then
        echo Source of genome build identified in the GTF file
    else
        echo Error: Source of genome build not identified in the GTF file
        exit 1;
    fi

    set -eo pipefail

    # Remove version suffix from transcript, gene, and exon IDs in order to match
    # previous Cell Ranger reference packages
    #
    # Input GTF:
    #     ... gene_id "ENSG00000223972.5"; ...
    # Output GTF:
    #     ... gene_id "ENSG00000223972"; gene_version "5"; ...

    gtf_modified="$(basename ~{annotation_gtf}).modified"
    # Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
    ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
    cat ~{annotation_gtf} \
        | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
        | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
        | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
        > "$gtf_modified"


    # Define string patterns for GTF tags
    # NOTES:
    # - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
    #   biotypes are part of a more generic "lncRNA" biotype.
    # - These filters are relevant only to GTF files from GENCODE. The GTFs from
    #   Ensembl release 98 have the following differences:
    #   - The names "gene_biotype" and "transcript_biotype" are used instead of
    #     "gene_type" and "transcript_type".
    #   - Readthrough transcripts are present but are not marked with the
    #     "readthrough_transcript" tag.
    #   - Only the X chromosome versions of genes in the pseudoautosomal regions
    #     are present, so there is no "PAR" tag.

    BIOTYPE_PATTERN=\
    "(protein_coding|lncRNA|\
    IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
    IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
    TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
    TR_V_pseudogene|TR_J_pseudogene)"
    GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
    TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
    READTHROUGH_PATTERN="tag \"readthrough_transcript\""
    PAR_PATTERN="tag \"PAR\""


    # Construct the gene ID allowlist. We filter the list of all transcripts
    # based on these criteria:
    #   - allowable gene_type (biotype)
    #   - allowable transcript_type (biotype)
    #   - no "PAR" tag (only present for Y chromosome PAR)
    #   - no "readthrough_transcript" tag
    # We then collect the list of gene IDs that have at least one associated
    # transcript passing the filters.

    cat "$gtf_modified" \
        | awk '$3 == "transcript"' \
        | grep -E "$GENE_PATTERN" \
        | grep -E "$TX_PATTERN" \
        | grep -Ev "$READTHROUGH_PATTERN" \
        | grep -Ev "$PAR_PATTERN" \
        | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
        | sort \
        | uniq \
        > "gene_allowlist"


    # Filter the GTF file based on the gene allowlist

    gtf_filtered="~{annotation_gtf_modified}"

    # Copy header lines beginning with "#"

    grep -E "^#" "$gtf_modified" > "$gtf_filtered"

    # Filter to the gene allowlist

    grep -Ff "gene_allowlist" "$gtf_modified" >> "$gtf_filtered"
    ls -lh *


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
    docker: "us.gcr.io/broad-gotc-prod/build-indices:1.0.0-2.7.10a-1683045573"
    memory: "50 GiB"
    disks: "local-disk ${disk} HDD"
    disk: disk + " GB" # TES
    cpu:"16"
  }
}

task BuildBWAreference {
  input {
    File genome_fa
    File? chrom_sizes_file

    # GTF annotation version refers to the version (GENCODE) or release (NCBI) listed in the GTF
    String gtf_annotation_version
    # Genome source can be NCBI or GENCODE
    String genome_source
    # Genome build is the assembly accession (NCBI) or version (GENCODE)
    String genome_build
    # Organism can be Macaque, Mouse, Human, etc.
    String organism
  }

String reference_name = "bwa-mem2-2.2.1-~{organism}-~{genome_source}-build-~{genome_build}"

  command <<<
    mkdir genome
    mv ~{chrom_sizes_file} genome/chrom.sizes
    file=~{genome_fa}
    if [ ${file: -3} == ".gz" ]
      then
      gunzip -c ~{genome_fa} > genome/genome.fa
    else
      mv ~{genome_fa} genome/genome.fa
    fi
    bwa-mem2 index genome/genome.fa
    tar --dereference -cvf - genome/ > ~{reference_name}.tar
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools-bwa-mem-2:1.0.0-2.2.1_x64-linux-1685469504"
    memory: "96GB"
    disks: "local-disk 100 HDD"
    disk: "100 GB" # TES
    cpu: "4"
  }

  output {
    File reference_bundle = "~{reference_name}.tar"
  }
}


