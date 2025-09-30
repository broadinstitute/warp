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

    Boolean run_add_introns = false
    Boolean run_mitofinder = false              
    String?  mito_accession                       # e.g. chimp or ferret mito accession (NC_…)
    File?    mito_ref_gbk                         # path to mitochondrion reference .gbk
    Array[String]? mitofinder_opts                # optional, override extra flags to MitoFinder/add_mito
  }

  # version of this pipeline
  String pipeline_version = "5.0.0"


  parameter_meta {
    annotations_gtf: "the annotation file"
    genome_fa: "the fasta file"
    biotypes: "gene_biotype attributes to include in the gtf file"
  }
  
    # ---- Append mitochondrial sequence + annotations ----
    # Note: String comparison is case-sensitive.
    if (run_mitofinder && organism != "Human" && organism != "Mouse") {
      call MitoAnnotate as mito {
        input:
          mito_accession = select_first([mito_accession]),
          mito_ref_gbk   = select_first([mito_ref_gbk]),
          genome_fa      = genome_fa,
          transcript_gtf = annotations_gtf,
          spec_name      = organism,
          mitofinder_opts = mitofinder_opts
      }
    }

    # Choose the files the rest of the pipeline should use:
    File genome_fa_for_indices = select_first([mito.out_fasta, genome_fa])
    File annotations_gtf_for_indices = select_first([mito.out_gtf, annotations_gtf])

    call BuildStarSingleNucleus {
      input:
        gtf_annotation_version = gtf_annotation_version,
        genome_fa = genome_fa_for_indices,
        annotation_gtf = annotations_gtf_for_indices,
        biotypes = biotypes,
        genome_build = genome_build,
        genome_source = genome_source,
        organism = organism
    }
    call CalculateChromosomeSizes {
      input:
        genome_fa = genome_fa_for_indices
    }
    call BuildBWAreference {
      input:
        genome_fa = genome_fa_for_indices,
        chrom_sizes_file = CalculateChromosomeSizes.chrom_sizes,
        genome_source = genome_source,
        genome_build = genome_build,
        gtf_annotation_version = gtf_annotation_version,
        organism = organism
    }

    call RecordMetadata {
      input:
        pipeline_version = pipeline_version,
        ### MODIFIED: Pass all mito-related info to the metadata task ###
        was_mitofinder_run = run_mitofinder,
        organism = organism,
        mito_accession_used = mito_accession,
        mito_ref_gbk_used = mito_ref_gbk,
        mitofinder_opts_used = mitofinder_opts,
        input_files = [annotations_gtf_for_indices, genome_fa_for_indices, biotypes],
        output_files = [
          BuildStarSingleNucleus.star_index,
          BuildStarSingleNucleus.modified_annotation_gtf,
          CalculateChromosomeSizes.chrom_sizes,
          BuildBWAreference.reference_bundle
        ]
    }

  if (run_add_introns) {
    call SNSS2AddIntronsToGTF {
      input:
      modified_annotation_gtf = BuildStarSingleNucleus.modified_annotation_gtf,
      genome_fa = genome_fa
    }
  }

  output {
    File snSS2_star_index = BuildStarSingleNucleus.star_index
    String pipeline_version_out = "BuildIndices_v~{pipeline_version}"
    File snSS2_annotation_gtf_modified = BuildStarSingleNucleus.modified_annotation_gtf
    File reference_bundle = BuildBWAreference.reference_bundle
    File chromosome_sizes = CalculateChromosomeSizes.chrom_sizes
    File metadata = RecordMetadata.metadata_file
    File? snSS2_annotation_gtf_with_introns = SNSS2AddIntronsToGTF.modified_annotation_gtf_with_introns
    # Optional outputs from the mito step
    File? mito_annotated_fasta = mito.out_fasta
    File? mito_annotated_gtf = mito.out_gtf
    File? star_index_with_introns = SNSS2AddIntronsToGTF.star_index_with_introns
  }
}

task MitoAnnotate {
  input {
    String mito_accession
    File   mito_ref_gbk
    File   genome_fa
    File   transcript_gtf
    String spec_name
    Array[String]? mitofinder_opts
  }

  command <<<
    set -euo pipefail

    # Run your Dockerized script; it writes *_mito.fasta and *_mito.gtf
    add_mito \
      -a ~{mito_accession} \
      -r ~{mito_ref_gbk} \
      -g ~{genome_fa} \
      -t ~{transcript_gtf} \
      -n ~{spec_name} \
      ~{sep=' ' mitofinder_opts}

    # List for debugging
    ls -lah
  >>>

  output {
    # add_mito names outputs by appending "_mito"
    # This regex-based approach is great and robust. No changes needed here.
    File out_fasta = sub(basename(genome_fa), "\\.(fa|fasta)(\\.gz)?$", "") + "_mito.fasta"
    File out_gtf   = sub(basename(transcript_gtf), "\\.(gtf)(\\.gz)?$", "") + "_mito.gtf"
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/add_mito:1.0.0-1.4.2"
    memory: "8 GiB"
    disks: "local-disk 50 HDD"
    cpu: 2
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
    # First check for marmoset GTF and modify header
    echo "checking for marmoset"
    if [[ "~{organism}" == "marmoset" || "~{organism}" == "Marmoset" ]]
    then
        echo "marmoset is detected, running header modification"
        python3 /script/create_marmoset_header_mt_genes.py \
            ~{annotation_gtf} > "/cromwell_root/header.gtf"
    else
        echo "marmoset is not detected"

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
    fi

    if [[ "~{organism}" == "marmoset" || "~{organism}" == "Marmoset" ]]
    then
        echo "marmoset detected, running marmoset GTF modification"
        echo "Listing files to check for head.gtf"
        ls
        python3 /script/modify_gtf_marmoset.py \
            --input-gtf "/cromwell_root/header.gtf" \
            --output-gtf ~{annotation_gtf_modified} \
            --species ~{organism}
        echo "listing files, should see modified gtf"
        ls
    else
        echo "running GTF modification for non-marmoset"
        python3 /script/modify_gtf.py \
            --input-gtf ~{annotation_gtf} \
            --output-gtf ~{annotation_gtf_modified} \
            --biotypes ~{biotypes}
    fi
    # python3 /script/modify_gtf.py  \
    # --input-gtf ~{annotation_gtf} \
    # --output-gtf ~{annotation_gtf_modified} \
    # --biotypes ~{biotypes}

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
    File modified_annotation_gtf = annotation_gtf_modified
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/build-indices:2.1.0"
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
    # New inputs for logging mito info
    Boolean was_mitofinder_run
    String organism
    String? mito_accession_used
    File? mito_ref_gbk_used
    Array[String]? mitofinder_opts_used
  }

  command <<<
    set -euo pipefail

    # create metadata file
    echo "Pipeline Version: ~{pipeline_version}" > metadata.txt
    echo "Date of Workflow Run: $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> metadata.txt
    echo "" >> metadata.txt

    ### NEW/MODIFIED: Add the MitoFinder section to the metadata file ###
    echo "--- MitoFinder Details ---" >> metadata.txt
    # Check if the boolean flag is true
    if [[ "~{was_mitofinder_run}" == "true" ]]; then
      echo "MitoFinder was run for organism: ~{organism}" >> metadata.txt
      # Check for and report the specific parameters used, handling optional inputs.
      if [ "~{mito_accession_used}" != "" ]; then
        echo "Mitochondrial Accession: ~{mito_accession_used}" >> metadata.txt
      fi
      if [ "~{mito_ref_gbk_used}" != "" ]; then
        echo "Mitochondrial Reference GBK: ~{mito_ref_gbk_used}" >> metadata.txt
      fi
      if [ "~{sep=' ' mitofinder_opts_used}" != "" ]; then
        echo "MitoFinder Extra Options: ~{sep=' ' mitofinder_opts_used}" >> metadata.txt
      fi
    else
      echo "MitoFinder was not run." >> metadata.txt
    fi
    echo "" >> metadata.txt
    ###################################################################

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

  task SNSS2AddIntronsToGTF {
  input {
    File modified_annotation_gtf
    File genome_fa
  }
    String basename = basename(modified_annotation_gtf, ".gtf")
    String star_index_name = "~{basename}_intron.tar"

  command <<<

  python3  /script/add-introns-to-gtf.py  \
    --input-gtf "~{modified_annotation_gtf}" \
    --output-gtf "~{basename}_with_introns.gtf"

  mkdir star
  STAR --runMode genomeGenerate \
  --genomeDir star \
  --genomeFastaFiles ~{genome_fa} \
  --sjdbGTFfile ~{basename}_with_introns.gtf \
  --sjdbOverhang 100 \
  --runThreadN 16

  tar -cvf ~{star_index_name} star
  >>>

  output {
    File modified_annotation_gtf_with_introns = "~{basename}_with_introns.gtf"
    File star_index_with_introns = star_index_name
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/build-indices:4.2.1"
    memory: "50 GiB"
    disks: "local-disk 100 HDD"
    disk: 100 + " GB" # TES
  }
}