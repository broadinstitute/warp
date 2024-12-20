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
  String pipeline_version = "3.2.0"


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

    call RecordMetadata {
      input:
      pipeline_version = pipeline_version,
      input_files = [annotations_gtf, genome_fa, biotypes],
      output_files = [
      BuildStarSingleNucleus.star_index,
      BuildStarSingleNucleus.modified_annotation_gtf,
      CalculateChromosomeSizes.chrom_sizes,
      BuildBWAreference.reference_bundle
      ]
    }

  output {
    File snSS2_star_index = BuildStarSingleNucleus.star_index
    String pipeline_version_out = "BuildIndices_v~{pipeline_version}"
    File snSS2_annotation_gtf_modified = BuildStarSingleNucleus.modified_annotation_gtf
    File reference_bundle = BuildBWAreference.reference_bundle
    File chromosome_sizes = CalculateChromosomeSizes.chrom_sizes
    File metadata = RecordMetadata.metadata_file
  }
}

task CalculateChromosomeSizes {
  input {
    File genome_fa
  }
  command <<<
    samtools faidx ~{genome_fa}
    cut -f1,2 "~{genome_fa}.fai" > chrom.sizes
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

  command <<<
    # Check that input GTF files contain input genome source, genome build version, and annotation version
    if head -10 ~{annotation_gtf} | grep -qi ~{genome_build}
    then
        echo "Genome version found in the GTF file"
    else
        echo "Error: Input genome version does not match version in GTF file"
        exit 1
    fi

    # Check that GTF file contains correct build source info in the first 10 lines of the GTF
    if head -10 ~{annotation_gtf} | grep -qi ~{genome_source}
    then
        echo "Source of genome build identified in the GTF file"
    else
        echo "Error: Source of genome build not identified in the GTF file"
        exit 1
    fi

    set -eo pipefail

    # Add gene_name attributes if missing
    if ! grep -q 'gene_name' ~{annotation_gtf}; then
        echo "Adding gene_name attributes to GTF"
        TEMP_GTF=$(mktemp)
        cat ~{annotation_gtf} | \
        awk -F'\t' 'BEGIN{OFS="\t"} {
            if ($1 ~ /^#/) {
                print
                next
            }
            attr=$9
            has_gene_name=0
            gene_val=""
            gene_id_val=""
            
            # Split attributes and check for gene_name
            split(attr, attrs, ";")
            for (i in attrs) {
                if (attrs[i] ~ /gene_name/) {
                    has_gene_name=1
                }
                if (attrs[i] ~ /gene /) {
                    match(attrs[i], /"([^"]+)"/)
                    gene_val=substr(attrs[i], RSTART+1, RLENGTH-2)
                }
                if (attrs[i] ~ /gene_id /) {
                    match(attrs[i], /"([^"]+)"/)
                    gene_id_val=substr(attrs[i], RSTART+1, RLENGTH-2)
                }
            }
            
            # If no gene_name, add it using gene or gene_id
            if (!has_gene_name) {
                if (gene_val != "") {
                    $9 = attr " gene_name \"" gene_val "\";"
                } else if (gene_id_val != "") {
                    $9 = attr " gene_name \"" gene_id_val "\";"
                }
            }
            print
        }' > "$TEMP_GTF"
        mv "$TEMP_GTF" ~{annotation_gtf}
    fi
    
    # Run modify_gtf.py first for all organisms
    python3 /script/modify_gtf.py \
    --input-gtf ~{annotation_gtf} \
    --output-gtf ~{annotation_gtf_modified} \
    --biotypes ~{biotypes}

    # Create a temporary final GTF path
    FINAL_GTF="~{annotation_gtf_modified}"
    
    # Conditionally process GTF for Marmoset
    if [[ "~{organism}" == "Marmoset" ]]; then
        echo "Processing Marmoset GTF file..."
        # Create a temporary GTF file path
        TEMP_GTF=$(mktemp)
        
        # Run GTF processing script  
        
        # Run the Marmoset-specific GTF processing on the modify_gtf.py output
        MARMOSET_OUTPUT=$(mktemp)
        bash /script/marmoset_gtf_.sh -i ~{annotation_gtf} -m ~{annotation_gtf_modified} -o "$MARMOSET_OUTPUT"        
        # Update the final GTF path
        mv "$MARMOSET_OUTPUT" ~{annotation_gtf_modified}
        FINAL_GTF="~{annotation_gtf_modified}"
    fi

    # Create STAR index using the appropriate GTF
    mkdir star
    STAR --runMode genomeGenerate \
    --genomeDir star \
    --genomeFastaFiles ~{genome_fa} \
    --sjdbGTFfile "$FINAL_GTF" \
    --sjdbOverhang 100 \
    --runThreadN 16

    tar -cvf ~{star_index_name} star
  >>>

  output {
    File star_index = star_index_name
    File modified_annotation_gtf = annotation_gtf_modified
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/build-indices:lk-fix-marmoset-buildindices"
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


task RecordMetadata {
  input {
    String pipeline_version
    Array[File] input_files
    Array[File] output_files
  }

  command <<<
    set -euo pipefail

    # create metadata file
    echo "Pipeline Version: ~{pipeline_version}" > metadata.txt
    echo "Date of Workflow Run: $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> metadata.txt
    echo "" >> metadata.txt

    # echo paths and md5sums for input files
    echo "Input Files and their md5sums:" >> metadata.txt
    for file in ~{sep=" " input_files}; do
      echo "$file : $(md5sum $file | awk '{print $1}')" >> metadata.txt
    done
    echo "" >> metadata.txt

    # echo paths and md5sums for input files
    echo "Output Files and their md5sums:" >> metadata.txt
    for file in ~{sep=" " output_files}; do
      echo "$file : $(md5sum $file | awk '{print $1}')" >> metadata.txt
    done
    echo "" >> metadata.txt

    # grab workspace bucket
    file="~{output_files[0]}"
    workspace_bucket=$(echo $file | awk -F'/' '{print $3}')
    echo "Workspace Bucket: $workspace_bucket" >> metadata.txt

    # grab submission ID
    submission_id=$(echo $file | awk -F'/' '{print $5}')
    echo "Submission ID: $submission_id" >> metadata.txt

    # grab workflow ID
    workflow_id=$(echo $file | awk -F'/' '{print $7}')
    echo "Workflow ID: $workflow_id" >> metadata.txt

    echo "" >> metadata.txt
  >>>

  output {
    File metadata_file = "metadata.txt"
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "5 GiB"
    disks: "local-disk 100 HDD"
    cpu: "1"
  }
}

