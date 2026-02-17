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
    File? biotypes                                # required only when run_modify_gtf=true for non-marmoset

    Boolean run_add_introns = false
    Boolean run_mitofinder = false
    Boolean run_modify_gtf

    String?  mito_accession                       # e.g. chimp or ferret mito accession (NC_…)
    File?    mito_ref_gbk                         # path to mitochondrion reference .gbk
    Array[String]? mitofinder_opts                # optional, override extra flags to MitoFinder/add_mito
    File?    annotations_gff                      # gff file for mitofinder
  }

  # version of this pipeline
  String pipeline_version = "5.1.0"

  parameter_meta {
    annotations_gtf: "the annotation file"
    genome_fa: "the fasta file"
    biotypes: "gene_biotype attributes to include in the gtf file; required when run_modify_gtf=true for non-marmoset organisms"
  }

  # ---- Mitofinder block (completely isolated; no effect when run_mitofinder = false) ----
  if (run_mitofinder) {
    call MitoAnnotate as annotate_with_mitofinder {
      input:
        mito_accession = select_first([mito_accession]), # check to make sure mito_accession is provided before running MitoAnnotate
        mito_ref_gbk   = select_first([mito_ref_gbk]), # check to make sure mito_ref_gbk is provided before running MitoAnnotate
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

    # Remove the old mito contig from the combined genome to avoid duplicates
    call RemoveDuplicateMitoContig {
      input:
        genome_fa      = annotate_with_mitofinder.out_fasta,
        mito_accession = select_first([mito_accession])
    }
  }

  # When run_mitofinder = true, use mito-processed files; otherwise pass through originals unchanged
  File final_genome_fa = select_first([RemoveDuplicateMitoContig.cleaned_fasta, genome_fa])
  File final_annotations_gtf = select_first([append_mito_gtf.out_gtf, annotations_gtf])

  # ---- Always-run GTF gene_name fix ----
  call FixGeneNames {
    input:
      annotation_gtf = final_annotations_gtf
  }

  # ---- Conditional GTF modification block ----
  Boolean is_marmoset = (organism == "marmoset" || organism == "Marmoset")

  # ---- Always-run GTF validation (non-marmoset only) ----
  # In the original monolithic BuildStarSingleNucleus, these checks ran
  # unconditionally for non-marmoset organisms—even when skip_gtf_modification
  # was true.  Extracting them here preserves that behaviour after the task
  # was split into separate steps.
  if (!is_marmoset) {
    call ValidateGTF {
      input:
        annotation_gtf = FixGeneNames.fixed_gtf,
        genome_source = genome_source,
        genome_build = genome_build
    }
  }

  if (run_modify_gtf && !is_marmoset) {
    call ModifyGTF {
      input:
        annotation_gtf = FixGeneNames.fixed_gtf,
        biotypes = select_first([biotypes])
    }
  }

  if (run_modify_gtf && is_marmoset) {
    call ModifyGTFMarmoset {
      input:
        annotation_gtf = FixGeneNames.fixed_gtf,
        organism = organism
    }
  }

  File gtf_for_star = select_first([ModifyGTF.modified_gtf, ModifyGTFMarmoset.modified_gtf, FixGeneNames.fixed_gtf])

  call BuildStarSingleNucleus {
    input:
      gtf_annotation_version = gtf_annotation_version,
      annotation_gtf = gtf_for_star,
      genome_fa = final_genome_fa,
      genome_build = genome_build,
      genome_source = genome_source,
      organism = organism,
      run_modify_gtf = run_modify_gtf
  }

  call CalculateChromosomeSizes {
    input:
      genome_fa = final_genome_fa
  }

  call BuildBWAreference {
    input:
      genome_fa = final_genome_fa,
      chrom_sizes_file = CalculateChromosomeSizes.chrom_sizes,
      genome_source = genome_source,
      genome_build = genome_build,
      gtf_annotation_version = gtf_annotation_version,
      organism = organism
  }

  # Centralize what goes into metadata
  # All optional outputs can be passed directly to select_all
  Array[File] recorded_inputs = select_all([
    annotations_gff,  # File? - already optional
    annotations_gtf,
    biotypes,
    genome_fa
  ])

  Array[File] recorded_outputs = select_all([
    annotate_with_mitofinder.out_fasta,  # File? from conditional block
    append_mito_gtf.out_gtf,              # File? from conditional block
    FixGeneNames.fixed_gtf,               # Always-run gene_name fix
    ModifyGTF.modified_gtf,               # File? from conditional block
    ModifyGTFMarmoset.modified_gtf,       # File? from conditional block
    BuildStarSingleNucleus.star_index,
    BuildStarSingleNucleus.modified_annotation_gtf,
    CalculateChromosomeSizes.chrom_sizes,
    BuildBWAreference.reference_bundle
  ])

  call RecordMetadata {
    input:
      pipeline_version = pipeline_version,
      organism = organism,
      genome_source = genome_source,
      genome_build = genome_build,
      gtf_annotation_version = gtf_annotation_version,
      run_mitofinder = run_mitofinder,
      run_modify_gtf = run_modify_gtf,
      run_add_introns = run_add_introns,
      is_marmoset = is_marmoset,
      input_files = recorded_inputs,
      output_files = recorded_outputs,
      input_annotations_gtf = annotations_gtf,
      input_genome_fa = genome_fa,
      mito_annotated_fasta = annotate_with_mitofinder.out_fasta,
      mito_appended_gtf = append_mito_gtf.out_gtf,
      fixed_gtf = FixGeneNames.fixed_gtf,
      modified_gtf = ModifyGTF.modified_gtf,
      modified_gtf_marmoset = ModifyGTFMarmoset.modified_gtf,
      star_annotation_gtf = BuildStarSingleNucleus.modified_annotation_gtf,
      star_index = BuildStarSingleNucleus.star_index
  }

  if (run_add_introns) {
    call SNSS2AddIntronsToGTF {
      input:
        modified_annotation_gtf = BuildStarSingleNucleus.modified_annotation_gtf,
        genome_fa = final_genome_fa
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

task RemoveDuplicateMitoContig {
  input {
    File genome_fa
    String mito_accession
  }

  command <<<
    set -euo pipefail

    cp ~{genome_fa} genome_input.fasta

    if grep -q "^>~{mito_accession}$" genome_input.fasta; then
      echo "Removing duplicate contig ~{mito_accession} from FASTA..."
      awk -v acc="~{mito_accession}" '
        BEGIN { deleted = 0 }
        $0 == ">" acc && deleted == 0 { deleted = 1; skip = 1; next }
        /^>/ { skip = 0 }
        !skip
      ' genome_input.fasta > cleaned_genome.fasta
    else
      echo "Contig ~{mito_accession} not found; no removal needed."
      cp genome_input.fasta cleaned_genome.fasta
    fi
  >>>

  output {
    File cleaned_fasta = "cleaned_genome.fasta"
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "4 GiB"
    disks: "local-disk 50 HDD"
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

task FixGeneNames {
  input {
    File annotation_gtf
  }

  meta {
    description: "Decompress GTF if needed and fix missing gene_name attributes by copying from gene_id"
  }

  command <<<
    set -eo pipefail

    # Decompress GTF if gzipped
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

    echo "GTF gene_name fix complete"
  >>>

  output {
    File fixed_gtf = "fixed_annotation.gtf"
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "4 GiB"
    disks: "local-disk 50 HDD"
    cpu: 1
  }
}

task ValidateGTF {
  input {
    File annotation_gtf
    String genome_source
    String genome_build
  }

  meta {
    description: "Validate that the GTF header contains the expected genome build and source. Runs for non-marmoset organisms regardless of run_modify_gtf."
  }

  command <<<
    set -eo pipefail

    GTF_FILE="~{annotation_gtf}"

    # Check that GTF contains expected genome build
    if head -10 ${GTF_FILE} | grep -qi ~{genome_build}
    then
        echo Genome version found in the GTF file
    else
        echo Error: Input genome version does not match version in GTF file
        exit 1;
    fi

    # Check that GTF contains expected genome source
    if head -10 ${GTF_FILE} | grep -qi ~{genome_source}
    then
        echo Source of genome build identified in the GTF file
    else
        echo Error: Source of genome build not identified in the GTF file
        exit 1;
    fi

    echo "GTF validation passed"
  >>>

  output {
    Boolean validation_passed = true
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "2 GiB"
    disks: "local-disk 10 HDD"
    cpu: 1
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
    Boolean run_modify_gtf
    Int disk = 100
  }

  meta {
    description: "Build reference index files for STAR aligner"
  }

  String gtf_prefix = if run_modify_gtf then "modified_" else ""
  String ref_name = "star2.7.10a-~{organism}-~{genome_source}-build-~{genome_build}-~{gtf_annotation_version}"
  String star_index_name = "~{gtf_prefix}~{ref_name}.tar"
  String annotation_gtf_modified = "~{gtf_prefix}v~{gtf_annotation_version}.annotation.gtf"

  command <<<
    # Decompress GTF if gzipped, otherwise copy to expected output name
    if [[ "~{annotation_gtf}" == *.gz ]]; then
        echo "Detected gzipped GTF file, decompressing..."
        gunzip -c ~{annotation_gtf} > ~{annotation_gtf_modified}
    else
        echo "GTF file is not compressed, copying..."
        cp ~{annotation_gtf} ~{annotation_gtf_modified}
    fi

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

task ModifyGTF {
  input {
    File annotation_gtf
    File biotypes
  }

  meta {
    description: "Modify GTF annotation file for non-marmoset organisms using biotype filtering"
  }

  command <<<
    set -eo pipefail

    GTF_FILE="~{annotation_gtf}"

    # Run standard GTF modification
    echo "Running GTF modification"
    python3 /script/modify_gtf.py \
        --input-gtf ${GTF_FILE} \
        --output-gtf modified.annotation.gtf \
        --biotypes ~{biotypes}
  >>>

  output {
    File modified_gtf = "modified.annotation.gtf"
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/build-indices:2.1.0"
    memory: "8 GiB"
    disks: "local-disk 50 HDD"
    cpu: 2
  }
}

task ModifyGTFMarmoset {
  input {
    File annotation_gtf
    String organism
  }

  meta {
    description: "Modify GTF annotation file for marmoset organisms"
  }

  command <<<
    set -eo pipefail

    GTF_FILE="~{annotation_gtf}"

    # Create marmoset header
    echo "Marmoset detected, running header modification"
    python3 /script/create_marmoset_header_mt_genes.py \
        ${GTF_FILE} > /cromwell_root/header.gtf

    # Run marmoset-specific GTF modification
    echo "Running marmoset GTF modification"
    python3 /script/modify_gtf_marmoset.py \
        --input-gtf /cromwell_root/header.gtf \
        --output-gtf modified.annotation.gtf \
        --species ~{organism}
  >>>

  output {
    File modified_gtf = "modified.annotation.gtf"
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/build-indices:2.1.0"
    memory: "8 GiB"
    disks: "local-disk 50 HDD"
    cpu: 2
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
    String organism
    String genome_source
    String genome_build
    String gtf_annotation_version
    Boolean run_mitofinder
    Boolean run_modify_gtf
    Boolean run_add_introns
    Boolean is_marmoset
    Array[File] input_files
    Array[File] output_files

    # Original inputs for tracking
    File input_annotations_gtf
    File input_genome_fa

    # Optional modification outputs for tracking which steps ran
    File? mito_annotated_fasta
    File? mito_appended_gtf
    File fixed_gtf
    File? modified_gtf
    File? modified_gtf_marmoset
    File star_annotation_gtf
    File star_index
  }

  command <<<
    set -euo pipefail

    # Helper: convert cromwell paths to gs:// paths
    to_gs() { echo "$1" | sed 's|^/mnt/disks/cromwell_root/|gs://|'; }

    # ---- Header ----
    echo "========================================" > metadata.txt
    echo "BuildIndices Pipeline Metadata" >> metadata.txt
    echo "========================================" >> metadata.txt
    echo "Pipeline Version: ~{pipeline_version}" >> metadata.txt
    echo "Date of Workflow Run: $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> metadata.txt
    echo "" >> metadata.txt

    # ---- Reference Genome Info ----
    echo "Reference Genome Configuration:" >> metadata.txt
    echo "  Organism: ~{organism}" >> metadata.txt
    echo "  Genome Source: ~{genome_source}" >> metadata.txt
    echo "  Genome Build: ~{genome_build}" >> metadata.txt
    echo "  GTF Annotation Version: ~{gtf_annotation_version}" >> metadata.txt
    echo "  Is Marmoset: ~{is_marmoset}" >> metadata.txt
    echo "" >> metadata.txt

    # ---- Pipeline Options ----
    echo "Pipeline Options:" >> metadata.txt
    echo "  run_mitofinder: ~{run_mitofinder}" >> metadata.txt
    echo "  run_modify_gtf: ~{run_modify_gtf}" >> metadata.txt
    echo "  run_add_introns: ~{run_add_introns}" >> metadata.txt
    echo "" >> metadata.txt

    # ---- Modifications Applied ----
    echo "Modifications Applied:" >> metadata.txt

    # MitoFinder
    if [ "~{run_mitofinder}" = "true" ]; then
      echo "  [MitoFinder] Ran mitochondrial annotation" >> metadata.txt
      echo "    Input genome FASTA: $(to_gs '~{input_genome_fa}')" >> metadata.txt
      if [ -n "~{default='NONE' mito_annotated_fasta}" ] && [ "~{default='NONE' mito_annotated_fasta}" != "NONE" ]; then
        echo "    Output mito-annotated FASTA: $(to_gs '~{mito_annotated_fasta}')" >> metadata.txt
      fi
      echo "    Input annotations GTF: $(to_gs '~{input_annotations_gtf}')" >> metadata.txt
      if [ -n "~{default='NONE' mito_appended_gtf}" ] && [ "~{default='NONE' mito_appended_gtf}" != "NONE" ]; then
        echo "    Output mito-appended GTF: $(to_gs '~{mito_appended_gtf}')" >> metadata.txt
      fi
    else
      echo "  [MitoFinder] Skipped" >> metadata.txt
    fi

    # FixGeneNames (always runs)
    echo "  [FixGeneNames] Fixed missing gene_name attributes" >> metadata.txt
    echo "    Output fixed GTF: $(to_gs '~{fixed_gtf}')" >> metadata.txt

    # GTF Modification
    if [ "~{run_modify_gtf}" = "true" ]; then
      if [ "~{is_marmoset}" = "true" ]; then
        echo "  [ModifyGTFMarmoset] Ran marmoset-specific GTF modification" >> metadata.txt
        echo "    Input GTF: $(to_gs '~{input_annotations_gtf}')" >> metadata.txt
        if [ -n "~{default='NONE' modified_gtf_marmoset}" ] && [ "~{default='NONE' modified_gtf_marmoset}" != "NONE" ]; then
          echo "    Output modified GTF: $(to_gs '~{modified_gtf_marmoset}')" >> metadata.txt
        fi
      else
        echo "  [ModifyGTF] Ran standard GTF modification" >> metadata.txt
        echo "    Input GTF: $(to_gs '~{input_annotations_gtf}')" >> metadata.txt
        if [ -n "~{default='NONE' modified_gtf}" ] && [ "~{default='NONE' modified_gtf}" != "NONE" ]; then
          echo "    Output modified GTF: $(to_gs '~{modified_gtf}')" >> metadata.txt
        fi
      fi
    else
      echo "  [ModifyGTF] Skipped" >> metadata.txt
    fi

    # STAR Index
    echo "  [BuildStarSingleNucleus] Built STAR index" >> metadata.txt
    echo "    Input GTF: $(to_gs '~{star_annotation_gtf}')" >> metadata.txt
    echo "    Output STAR index: $(to_gs '~{star_index}')" >> metadata.txt

    echo "" >> metadata.txt

    # ---- Input Files ----
    echo "Input Files and their md5sums:" >> metadata.txt
    for file in ~{sep=" " input_files}; do
      gs_path=$(to_gs "$file")
      echo "  $gs_path : $(md5sum "$file" | awk '{print $1}')" >> metadata.txt
    done
    echo "" >> metadata.txt

    # ---- Output Files ----
    echo "Output Files and their md5sums:" >> metadata.txt
    for file in ~{sep=" " output_files}; do
      gs_path=$(to_gs "$file")
      echo "  $gs_path : $(md5sum "$file" | awk '{print $1}')" >> metadata.txt
    done
    echo "" >> metadata.txt

    # ---- Cromwell Execution Info ----
    file="~{output_files[0]}"
    workspace_bucket=$(echo $file | awk -F'/' '{print $3}')
    echo "Workspace Bucket: $workspace_bucket" >> metadata.txt

    submission_id=$(echo $file | awk -F'/' '{print $5}')
    echo "Submission ID: $submission_id" >> metadata.txt

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
