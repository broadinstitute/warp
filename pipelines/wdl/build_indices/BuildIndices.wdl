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
    Boolean skip_gtf_modification = false

    String?  mito_accession                       # e.g. chimp or ferret mito accession (NC_…)
    File?    mito_ref_gbk                         # path to mitochondrion reference .gbk
    Array[String]? mitofinder_opts                # optional, override extra flags to MitoFinder/add_mito
    File?    annotations_gff                       # gff file for mitofinder
  }

  if (false) {
    String? none = "None"
  }

  # version of this pipeline
  String pipeline_version = "5.0.3"


  parameter_meta {
    annotations_gtf: "the annotation file"
    genome_fa: "the fasta file"
    biotypes: "gene_biotype attributes to include in the gtf file"
  }

    # ---- Append mitochondrial sequence + annotations ----
    # Note: String comparison is case-sensitive.
    if (run_mitofinder) {
      call MitoAnnotate as annotate_with_mitofinder {
        input:
          mito_accession = select_first([mito_accession]),
          mito_ref_gbk   = select_first([mito_ref_gbk]),
          genome_fa      = genome_fa,
          transcript_gff = select_first([annotations_gff]),
          spec_name      = organism,
          mitofinder_opts = mitofinder_opts
      }
      call AppendMitoGTF as append_mito_gtf {
        input:
          original_gtf   = annotations_gtf,
          mito_gtf       = annotate_with_mitofinder.out_gtf
      }
    }

    call BuildStarSingleNucleus {
      input:
        gtf_annotation_version = gtf_annotation_version,
        annotation_gtf = select_first([if run_mitofinder then append_mito_gtf.out_gtf else none, annotations_gtf]),
        genome_fa = genome_fa, #select_first([if run_mitofinder then annotate_with_mitofinder.out_fasta else none, genome_fa]),
        biotypes = biotypes,
        genome_build = genome_build,
        genome_source = genome_source,
        organism = organism,
        skip_gtf_modification = skip_gtf_modification,
        mito_accession = select_first([mito_accession])
    }
    call CalculateChromosomeSizes {
      input:
        genome_fa = select_first([if run_mitofinder then annotate_with_mitofinder.out_fasta else none, genome_fa]),
    }
    call BuildBWAreference {
      input:
        genome_fa = select_first([if run_mitofinder then annotate_with_mitofinder.out_fasta else none, genome_fa]),
        chrom_sizes_file = CalculateChromosomeSizes.chrom_sizes,
        genome_source = genome_source,
        genome_build = genome_build,
        gtf_annotation_version = gtf_annotation_version,
        organism = organism,
        mito_accession = select_first([mito_accession])
    }

    call RecordMetadata {
      input:
        pipeline_version = pipeline_version,
        was_mitofinder_run = run_mitofinder,
        organism = organism,
        mito_accession_used = mito_accession,
        mito_ref_gbk_used = mito_ref_gbk,
        mitofinder_opts_used = mitofinder_opts,
        input_files = select_all([if run_mitofinder then annotations_gff else none, annotations_gtf, biotypes, genome_fa]),
        output_files = select_all([
                                  if run_mitofinder then annotate_with_mitofinder.out_fasta else none,
                                  if run_mitofinder then append_mito_gtf.out_gtf else none,
                                  BuildStarSingleNucleus.star_index,
                                  BuildStarSingleNucleus.modified_annotation_gtf,
                                  CalculateChromosomeSizes.chrom_sizes,
                                  BuildBWAreference.reference_bundle
                                  ])
    }

  if (run_add_introns) {
    call SNSS2AddIntronsToGTF {
      input:
      modified_annotation_gtf = BuildStarSingleNucleus.modified_annotation_gtf,
      genome_fa = genome_fa
    }
  }

  output {
    File star_index_tar = BuildStarSingleNucleus.star_index
    String pipeline_version_out = "BuildIndices_v~{pipeline_version}"
    File snSS2_annotation_gtf_modified = BuildStarSingleNucleus.modified_annotation_gtf
    File reference_bundle = BuildBWAreference.reference_bundle
    File chromosome_sizes = CalculateChromosomeSizes.chrom_sizes
    File metadata = RecordMetadata.metadata_file
    File? snSS2_annotation_gtf_with_introns = SNSS2AddIntronsToGTF.modified_annotation_gtf_with_introns
    # Optional outputs from the mito step
    File? mito_annotated_fasta = annotate_with_mitofinder.out_fasta
    File? mito_annotated_gtf = append_mito_gtf.out_gtf
    File? star_index_with_introns = SNSS2AddIntronsToGTF.star_index_with_introns
  }
}

task MitoAnnotate {
  input {
    String mito_accession
    File   mito_ref_gbk
    File   genome_fa
    File   transcript_gff
    String spec_name
    Array[String]? mitofinder_opts
  }

  command <<<
    set -euo pipefail
    set -x

    # WDL/Cromwell already localizes File inputs; these are valid inside the container.
    GBK="~{mito_ref_gbk}"
    FA="~{genome_fa}"
    GFF_IN="~{transcript_gff}"

    # Make a temp dir in the container for any conversions/decompression
    tmpDir="$(mktemp -d)"
    chmod 777 "$tmpDir"

    echo "Container PWD: $(pwd)"
    echo "Checking mounted /cromwell_root (should exist in Cromwell docker backends):"
    ls -ld /cromwell_root || true

    echo "Preflight: show inputs"
    ls -l "$GBK" || (echo "Missing GBK at $GBK" && exit 1)
    ls -l "$FA"  || (echo "Missing FASTA at $FA" && exit 1)
    ls -l "$GFF_IN" || (echo "Missing transcript file at $GFF_IN" && exit 1)

    GFF="$GFF_IN"

    # If add_mito can't read .gz,  decompress to a working copy
    if gzip -t "$GFF_IN" 2>/dev/null; then
      echo "Decompressing transcript annotations to GFF (stream-safe)..."
      zcat -f "$GFF_IN" > "$tmpDir/transcripts.gff"
      GFF="$tmpDir/transcripts.gff"
    fi


    echo "Running add_mito with:"

    # Run add_mito using ONLY the localized, in-container paths
    add_mito \
      -a "~{mito_accession}" \
      -r "$GBK" \
      -g "$FA" \
      -t "$GFF" \
      -n "~{spec_name}" \
      ~{sep=' ' mitofinder_opts}


    FA_OUT=$(find /mnt/disks/cromwell_root -name '*_mito.fasta' | head -n1)
    GTF_OUT=$(find /mnt/disks/cromwell_root -name '*_mito.gtf' | head -n1)

    cp "$FA_OUT" ./genome_mito.fasta
    cp "$GTF_OUT" ./transcripts_mito.gtf

  >>>

  output {
    File out_fasta = "genome_mito.fasta"
    File out_gtf   = "transcripts_mito.gtf"
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/add_mito:1.0.0-1.4.2"
    memory: "8 GiB"
    disks: "local-disk 50 HDD"
    cpu: 2
  }
}

task AppendMitoGTF {
  input {
    File original_gtf
    File mito_gtf
  }

  command <<<
    set -euo pipefail

    # grep mitofinder in the mito_gtf and append those lines to the original gtf
    grep "mitofinder" ~{mito_gtf} > mito_only.gtf

    # Check if original GTF is gzipped and decompress if needed
    if [[ "~{original_gtf}" == *.gz ]]; then
      echo "Original GTF is gzipped, decompressing..."
      gunzip -c ~{original_gtf} > original_decompressed.gtf
      ORIGINAL_FILE="original_decompressed.gtf"
    else
      ORIGINAL_FILE="~{original_gtf}"
    fi

    # Concatenate the original GTF and the mito GTF
    echo "Combining GTF files..."
    cat "${ORIGINAL_FILE}" mito_only.gtf > combined_annotations.gtf

  >>>

  output {
    File out_gtf = "combined_annotations.gtf"
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "2 GiB"
    disks: "local-disk 10 HDD"
    cpu: 1
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
    Boolean skip_gtf_modification
    Int disk = 100
    String? mito_accession
  }

  meta {
    description: "Modify GTF files and build reference index files for STAR aligner"
  }

  String gtf_prefix = if skip_gtf_modification then "" else "modified_"
  String ref_name = "star2.7.10a-~{organism}-~{genome_source}-build-~{genome_build}-~{gtf_annotation_version}"
  String star_index_name = "~{gtf_prefix}~{ref_name}.tar"
  String annotation_gtf_modified = "~{gtf_prefix}v~{gtf_annotation_version}.annotation.gtf"

  command <<<
    # Decompress GTF if it's gzipped
    if [[ "~{annotation_gtf}" == *.gz ]]; then
        echo "Detected gzipped GTF file, decompressing..."
        gunzip -c ~{annotation_gtf} > annotation.gtf
        GTF_FILE="annotation.gtf"
    else
        echo "GTF file is not compressed"
        GTF_FILE="~{annotation_gtf}"
    fi

    # Fix missing gene_name attributes
    echo "Checking and fixing gene_name attributes in GTF..."
    awk -F'\t' 'BEGIN { OFS="\t" }
      /^#/ { print; next }
      {
        gene_id = ""; gene_name = "";
        if ($9 ~ /gene_id/) {
          n = split($9, a, /gene_id "/)
          if (n > 1) {
            split(a[2], b, "\"")
            gene_id = b[1]
          }
        }

        # Check if gene_name is missing and add it
        if ($9 !~ /gene_name/ && gene_id != "") {
          sub(/[[:space:]]*;[[:space:]]*$/, "", $9)  # remove trailing semicolons/spaces
          $9 = $9 "; gene_name \"" gene_id "\";"
        }

        print
      }' "$GTF_FILE" > fixed_annotation.gtf

    # Use the fixed GTF for downstream processing
    GTF_FILE="fixed_annotation.gtf"
    echo "GTF gene_name fix complete"

    # First check for marmoset GTF and modify header
    echo "checking for marmoset"
    if [[ "~{organism}" == "marmoset" || "~{organism}" == "Marmoset" ]]
    then
        echo "marmoset is detected, running header modification"
        python3 /script/create_marmoset_header_mt_genes.py \
            ${GTF_FILE} > "/cromwell_root/header.gtf"
    else
        echo "marmoset is not detected"

        # Check that input GTF files contain input genome source, genome build version, and annotation version
        if head -10 ${GTF_FILE} | grep -qi ~{genome_build}
        then
            echo Genome version found in the GTF file
        else
            echo Error: Input genome version does not match version in GTF file
            exit 1;
        fi

        # Check that GTF file contains correct build source info in the first 10 lines of the GTF
        if head -10 ${GTF_FILE} | grep -qi ~{genome_source}
        then
            echo Source of genome build identified in the GTF file
        else
            echo Error: Source of genome build not identified in the GTF file
            exit 1;
        fi
        set -eo pipefail
    fi

    if [ "~{skip_gtf_modification}" = "false" ]; then
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
                --input-gtf ${GTF_FILE} \
                --output-gtf ~{annotation_gtf_modified} \
                --biotypes ~{biotypes}
        fi
    else
        echo "Skipping GTF modification — using original GTF for STAR index"
        cp ${GTF_FILE} ~{annotation_gtf_modified}
    fi

    ## --- Remove duplicate mito contig if mito_accession is set
    #if [ -n "~{mito_accession}" ]; then
    #  echo "mito_accession provided: ~{mito_accession}"
#
    #  if grep -q "^>~{mito_accession}$" ~{genome_fa}; then
    #    echo "Removing duplicate contig ~{mito_accession} from FASTA..."
#
    #    awk -v acc="~{mito_accession}" '
    #      BEGIN { deleted = 0 }
    #      $0 == ">" acc && deleted == 0 { deleted = 1; skip = 1; next }
    #      /^>/ { skip = 0 }
    #      !skip
    #      ' ~{genome_fa} > genome_mito.filtered.fasta
    #    mv genome_mito.filtered.fasta ~{genome_fa}
    #  else
    #    echo "Contig ~{mito_accession} not found; skipping removal."
    #  fi
    #else
    #    echo "No mito_accession provided, skipping contig removal."
    #fi

    mkdir star
    STAR --runMode genomeGenerate \
    --genomeDir star \
    --genomeFastaFiles ~{genome_fa} \
    --sjdbGTFfile ~{annotation_gtf_modified} \
    --sjdbOverhang 100 \
    --runThreadN 16 \
    --limitGenomeGenerateRAM=43375752629

    tar -cvf ~{star_index_name} star
  >>>

  output {
    File star_index = star_index_name
    File modified_annotation_gtf = annotation_gtf_modified
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/build-indices:2.1.0"
    memory: "64 GiB"
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
    String? mito_accession
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

    # --- Remove duplicate contig if mito_accession is provided ---
    if [ -n "~{mito_accession}" ]; then
      echo "mito_accession provided: ~{mito_accession}"

      if grep -q "^>~{mito_accession}$" genome/genome.fa; then
        echo "Removing duplicate contig ~{mito_accession} from FASTA..."

        awk -v acc="~{mito_accession}" '
          BEGIN { deleted = 0 }
          $0 == ">" acc && deleted == 0 { deleted = 1; skip = 1; next }
          /^>/ { skip = 0 }
          !skip
        ' genome/genome.fa > genome/genome.filtered.fa

        mv genome/genome.filtered.fa genome/genome.fa
      else
        echo "Contig ~{mito_accession} not found in genome.fa, no removal needed."
      fi
    else
        echo "No mito_accession provided, skipping contig removal."
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

    echo "--- MitoFinder Details ---" >> metadata.txt
    # Check if the boolean flag is true
    if [[ "~{was_mitofinder_run}" == "true" ]]; then
      echo "MitoFinder was run for organism: ~{organism}" >> metadata.txt
      # Check for and report the specific parameters used, handling optional inputs.
      if [ "~{mito_accession_used}" != "" ]; then
        echo "Mitochondrial Accession: ~{mito_accession_used}" >> metadata.txt
      fi
      if [ "~{mito_ref_gbk_used}" != "" ]; then
        echo "Mitochondrial Reference GBK: ~{mito_ref_gbk_used} (md5sum: $(md5sum "~{mito_ref_gbk_used}" | awk '{print $1}'))" >> metadata.txt
      fi
      if [ "~{sep=' ' mitofinder_opts_used}" != "" ]; then
        echo "MitoFinder Extra Options: ~{sep=' ' mitofinder_opts_used}" >> metadata.txt
      fi
    else
      echo "MitoFinder was not run." >> metadata.txt
    fi
    echo "" >> metadata.txt

    # echo paths and md5sums for input files
    echo "Input Files and their md5sums:" >> metadata.txt
    for file in ~{sep=" " input_files}; do
      gs_path=$(echo "$file" | sed 's|^/mnt/disks/cromwell_root/|gs://|')
      echo "$gs_path : $(md5sum "$file" | awk '{print $1}')" >> metadata.txt
    done
    echo "" >> metadata.txt

    # echo paths and md5sums for input files
    echo "Output Files and their md5sums:" >> metadata.txt
    for file in ~{sep=" " output_files}; do
      gs_path=$(echo "$file" | sed 's|^/mnt/disks/cromwell_root/|gs://|')
      echo "$gs_path : $(md5sum "$file" | awk '{print $1}')" >> metadata.txt
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
