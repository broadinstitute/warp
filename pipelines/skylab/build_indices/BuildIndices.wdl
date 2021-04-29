version 1.0

struct References {
  File genome_fa
  File annotation_gtf
}

task GetReferences {
  input {
    String gtf_version
    String organism
    String organism_prefix
  }

  meta {
    description: "Download files needed for building the designated references"
  }

  String ftp_path = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_~{organism}/release_~{gtf_version}"
  String genome_fa = "GRC~{organism_prefix}38.primary_assembly.genome.fa"
  String annotation_gtf = "gencode.v~{gtf_version}.primary_assembly.annotation.gtf"

  command <<<
    set -eo pipefail

    ## download fasta
    wget ~{ftp_path}/~{genome_fa}.gz
    gunzip ~{genome_fa}.gz

    ## download gtf file
    wget ~{ftp_path}/~{annotation_gtf}.gz
    gunzip ~{annotation_gtf}.gz
  >>>

  output {
      References references = object {
      genome_fa: genome_fa,
      annotation_gtf: annotation_gtf
    }
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-star:v2.7.8a"
    disks: "local-disk 10 HDD"
  }
}

task BuildStar {
  input {
    String gtf_version
    String organism
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
    docker: "quay.io/humancellatlas/secondary-analysis-star:v2.7.8a"
    memory: "50 GiB"
    disks :"local-disk 100 HDD"
    cpu:"16"
  }
}

task BuildStarSingleNucleus {
  input {
    String gtf_version
    String organism
    References references
  }

  meta {
    description: "Modify gtf files and build reference index files for STAR aligner"
  }

  String ref_name = "star_primary_gencode_~{organism}_v~{gtf_version}"
  String star_index_name = "~{ref_name}_modified.tar"
  String genome_fa = "star_primary_gencode_~{organism}_v~{gtf_version}"
  String genome_fa_modified = "modified_star_primary_gencode_~{organism}_v~{gtf_version}"
  String annotation_gtf_modified = "modified_gencode.v~{gtf_version}.primary_assembly.annotation.gtf"
  command <<<
    set -eo pipefail
    if ~{organism} == "mouse"
      then
        # sed commands:
        # 1. Replace metadata after space with original contig name, as in GENCODE
        # 2. Add "chr" to names of autosomes and sex chromosomes
        # 3. Handle the mitochrondrial chromosome
        cat ~{genome_fa} \
        | sed -E 's/^>(\S+).*/>\1 \1/' \
        | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
        | sed -E 's/^>MT />chrM /' \
        > genome_fa_modfied_file


        # Remove version suffix from transcript, gene, and exon IDs in order to match
        # previous Cell Ranger reference packages
        #
        # Input GTF:
        #     ... gene_id "ENSMUSG00000102693.1"; ...
        # Output GTF:
        #     ... gene_id "ENSMUSG00000102693"; gene_version "1"; ...
        gtf_modified="$build/$(basename "$gtf_in").modified"
        # Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
        ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
        cat "$gtf_in" \
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
        BIOTYPE_PATTERN=\
        "(protein_coding|lncRNA|\
        IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
        IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
        TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
        TR_V_pseudogene|TR_J_pseudogene)"
        GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
        TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
        READTHROUGH_PATTERN="tag \"readthrough_transcript\""


        # Construct the gene ID allowlist. We filter the list of all transcripts
        # based on these criteria:
        #   - allowable gene_type (biotype)
        #   - allowable transcript_type (biotype)
        #   - no "readthrough_transcript" tag
        # We then collect the list of gene IDs that have at least one associated
        # transcript passing the filters.
        cat "$gtf_modified" \
        | awk '$3 == "transcript"' \
        | grep -E "$GENE_PATTERN" \
        | grep -E "$TX_PATTERN" \
        | grep -Ev "$READTHROUGH_PATTERN" \
        | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
        | sort \
        | uniq \
        > "${build}/gene_allowlist"


        # Filter the GTF file based on the gene allowlist
        gtf_filtered="${build}/$(basename "$gtf_in").filtered"
        # Copy header lines beginning with "#"
        grep -E "^#" "$gtf_modified" > "$gtf_filtered"
        # Filter to the gene allowlist
        grep -Ff "${build}/gene_allowlist" "$gtf_modified" \
        >> "$gtf_filtered"
    fi
    if ~{organism} == "human"
    then
    # Genome metadata
    genome="GRCh38"
    version="2020-A"


    # Set up source and build directories
    build="GRCh38-2020-A_build"
    mkdir -p "$build"


    # Download source files if they do not exist in reference_sources/ folder
    source="reference_sources"
    mkdir -p "$source"


    fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    fasta_in="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"
    gtf_in="${source}/gencode.v32.primary_assembly.annotation.gtf"


    if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" | zcat > "$fasta_in"
    fi
    if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" | zcat > "$gtf_in"
    fi


    # Modify sequence headers in the Ensembl FASTA to match the file
    # "GRCh38.primary_assembly.genome.fa" from GENCODE. Unplaced and unlocalized
    # sequences such as "KI270728.1" have the same names in both versions.
    #
    # Input FASTA:
    #   >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
    #
    # Output FASTA:
    #   >chr1 1
    fasta_modified="$build/$(basename "$fasta_in").modified"
    # sed commands:
    # 1. Replace metadata after space with original contig name, as in GENCODE
    # 2. Add "chr" to names of autosomes and sex chromosomes
    # 3. Handle the mitochrondrial chromosome
    cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"


    # Remove version suffix from transcript, gene, and exon IDs in order to match
    # previous Cell Ranger reference packages
    #
    # Input GTF:
    #     ... gene_id "ENSG00000223972.5"; ...
    # Output GTF:
    #     ... gene_id "ENSG00000223972"; gene_version "5"; ...
    gtf_modified="$build/$(basename "$gtf_in").modified"
    # Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
    ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
    cat "$gtf_in" \
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
    > "${build}/gene_allowlist"


    # Filter the GTF file based on the gene allowlist
    gtf_filtered="${build}/$(basename "$gtf_in").filtered"
    # Copy header lines beginning with "#"
    grep -E "^#" "$gtf_modified" > "$gtf_filtered"
    # Filter to the gene allowlist
    grep -Ff "${build}/gene_allowlist" "$gtf_modified" \
    >> "$gtf_filtered"

    fi

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
    File genome_fa_modfied = genome_fa_modfied_file
    References references = object {
             genome_fa: genome_fa_modified,
           annotation_gtf: annotation_gtf_modified
           }
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-star:v2.7.8a"
    memory: "50 GiB"
    disks :"local-disk 100 HDD"
    cpu:"16"
  }
}

task BuildRsem {
  input {
    String gtf_version
    String organism
    References references
  }

  meta {
    description: "build reference index files for RSEM"
  }

  String ref_name = "rsem_primary_gencode_~{organism}_v~{gtf_version}"
  String rsem_index_name = "~{ref_name}.tar"

  command {
    set -eo pipefail
    mkdir rsem
    rsem-prepare-reference --gtf ~{references.annotation_gtf} --bowtie ~{references.genome_fa} rsem/rsem_trans_index
    tar -cvf ~{rsem_index_name} rsem/
  }
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-rsem:v0.2.2-1.3.0"
    memory: "10 GiB"
    disks: "local-disk 100 HDD"
  }
  output {
    File rsem_index = rsem_index_name
  }
}

task BuildHisat2FromRsem {
  input {
    File rsem_index
  }

  String rsem_reference_name = basename(rsem_index, ".tar")
  String ref_name = "hisat2_from_rsem_~{rsem_reference_name}"
  String hisat2_index_name = "~{ref_name}.tar.gz"

  command {

    # extract rsem index
    tar -xf ~{rsem_index}

    # build index
    hisat2-build -p 8 rsem/rsem_trans_index.idx.fa ~{ref_name}
    mkdir ~{ref_name}
    mv ./*.ht2 ~{ref_name}
    tar -zcvf ~{hisat2_index_name} ~{ref_name}
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory: "8 GiB"
    disks: "local-disk 100 HDD"
    cpu: "8"
  }

  output {
    File hisat2_index = hisat2_index_name
  }
}

task BuildHisat2 {
  # This method builds a HISAT2 index directly from the downloaded GTF and Fasta files on GENCODE
  # It does not include dbsnp, whereas the below method does.

  input {
    String gtf_version
    String organism
    File genome_fa
  }

  String ref_name = "hisat2_primary_gencode_~{organism}_v~{gtf_version}"
  String hisat2_index_name = "~{ref_name}.tar.gz"

  command {
    set -eo pipefail

    # build index
    hisat2-build -p 8 ~{genome_fa} ~{ref_name}
    mkdir ~{ref_name}
    mv ./*.ht2 ~{ref_name}
    tar -zcvf ~{hisat2_index_name} ~{ref_name}
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory: "64 GiB"
    disks: "local-disk 100 HDD"
    cpu: "8"
  }

  output {
    File hisat2_index = hisat2_index_name
  }
}

task BuildHisat2SnpHaplotypeSplicing {
  # This version includes SNP, haplotype, and splices
  input {
    String organism
    String genome_short_string
    References references

    String gtf_version
    String dbsnp_version
  }

  String ref_name = "star_primary_gencode_~{organism}_v~{gtf_version}"
  String snp_file = "snp~{dbsnp_version}Common.txt"
  String hisat2_index_name = "~{ref_name}.tar.gz"

  command <<<

    HISAT2_DIR=/opt/tools/hisat2-2.1.0

    # Compressed fasta required here
    gzip ~{references.genome_fa}

    # download snp file
    wget http://hgdownload.cse.ucsc.edu/goldenPath/~{genome_short_string}/database/~{snp_file}.gz
    gunzip ~{snp_file}.gz

    # extract snps, splice sites, and exon information
    $HISAT2_DIR/hisat2_extract_snps_UCSC.py ~{references.genome_fa}.gz ~{snp_file} genome
    $HISAT2_DIR/hisat2_extract_splice_sites.py ~{references.annotation_gtf} > genome.ss
    $HISAT2_DIR/hisat2_extract_exons.py ~{references.annotation_gtf} > genome.exon

    # build the hisat2 reference
    $HISAT2_DIR/hisat2-build \
      -p 8 \
      genome.fa \
      --snp genome.snp \
      --haplotype genome.haplotype \
      --ss genome.ss \
      --exon genome.exon \
      genome_snp_tran

    mkdir ~{ref_name}
    cp ./*.ht2 ~{ref_name}
    tar -zcvf ~{hisat2_index_name} ~{ref_name}

  >>>

  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.3.0-2-2.1.0"
    memory: "240 GiB"
    disks: "local-disk 100 HDD"
    cpu: "16"
  }
  output {
    File hisat2_index = hisat2_index_name
  }
}

task BuildPicardRefFlat {
  input {
    References references
  }

  String refflat_name = basename(references.annotation_gtf, ".gtf") + ".refflat.txt"

  command {
    set -eo pipefail

    gtfToGenePred -genePredExt -geneNameAsName2  ~{references.annotation_gtf} refflat.tmp.txt

    paste <(cut -f 12 refflat.tmp.txt) <(cut -f 1-10 refflat.tmp.txt) > ~{refflat_name}

  }

  runtime {
    docker: "quay.io/humancellatlas/gtf_to_genepred:v0.0.0"
    memory: "8 GiB"
    disks: "local-disk 100 HDD"
    cpu: "8"
  }

  output {
      File refflat = refflat_name
  }
}

task BuildIntervalList {
  input {
    References references
  }

  String interval_list_name = basename(references.annotation_gtf, ".gtf") + ".interval_list"

  command <<<
    set -eo pipefail


    # index the fasta file

    samtools faidx ~{references.genome_fa}
    cut -f1,2 ~{references.genome_fa}.fai > sizes.genome

    awk -F '\t'  '{  printf "@SQ\tSN:%s\tLN:%s\n", $1, $2 }' sizes.genome  >> ~{interval_list_name}

    grep 'gene_type "rRNA"' ~{references.annotation_gtf} |
        awk '$3 == "transcript"' |
    cut -f1,4,5,7,9 |
    perl -lane '
        /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
        print join "\t", (@F[0,1,2,3], $1)
    ' |
    sort -k1V -k2n -k3n  >> ~{interval_list_name}

  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-umitools:0.0.1"
    memory: "8 GiB"
    disks: "local-disk 100 HDD"
    cpu: "8"
  }

  output {
      File interval_list = interval_list_name
  }
}


workflow BuildIndices {
  input {
    String gtf_version
    String organism
    String organism_prefix
    String genome_short_string
    String dbsnp_version
  }

  # version of this pipeline

  String pipeline_version = "0.0.1"

  parameter_meta {
    gtf_version: "the actual number of gencode, ex.  27"
    organism: "Either 'human' or 'mouse'"
    organism_prefix: "Either 'h' or 'm'"
    genome_short_string: "e.g. hg38, mm10"
    dbsnp_version: "integer num, ex 150"
  }

  call GetReferences {
    input:
      gtf_version = gtf_version,
      organism = organism,
      organism_prefix = organism_prefix
  }

  call BuildPicardRefFlat   {
     input:
        references = GetReferences.references
  }

  call BuildIntervalList  {
      input:
        references = GetReferences.references
  }

  call BuildStar {
    input:
      gtf_version = gtf_version,
      organism = organism,
      references = GetReferences.references
  }

  call BuildRsem {
    input:
      gtf_version = gtf_version,
      organism = organism,
      references = GetReferences.references
  }

  call BuildHisat2FromRsem {
    input:
      rsem_index = BuildRsem.rsem_index
  }

  call BuildHisat2 {
    input:
      gtf_version = gtf_version,
      organism = organism,
      genome_fa = GetReferences.references.genome_fa
  }

  call BuildHisat2SnpHaplotypeSplicing {
    input:
      gtf_version = gtf_version,
      organism = organism,
      genome_short_string = genome_short_string,
      dbsnp_version = dbsnp_version,
      references = GetReferences.references
  }

  output {
    File star_index = BuildStar.star_index
    File rsem_index = BuildRsem.rsem_index
    File hisat2_from_rsem_index = BuildHisat2FromRsem.hisat2_index
    File hisat2_index = BuildHisat2.hisat2_index
    File hisat2_snp_haplotype_splicing_index = BuildHisat2SnpHaplotypeSplicing.hisat2_index
    File refflat = BuildPicardRefFlat.refflat
    File interval_list = BuildIntervalList.interval_list

    File genome_fa = GetReferences.references.genome_fa
    File annotation_gtf = GetReferences.references.annotation_gtf
  }
}
