version 1.0

workflow PrepareTransWindow {
  input {
    String window_id
    File genome_dosage
    File genome_dosage_tbi
    File trans_window_associations
    File expression_phenotypes
    File splicing_phenotypes
    File protein_phenotypes
    File target_phenotypes
    File? expression_phenotypes_tbi
    File? expression_phenotype_lookup
    File? splicing_phenotypes_tbi
    File? splicing_phenotype_lookup
    File? protein_phenotypes_tbi
    File? protein_phenotype_lookup
    Int top_n_expression = 25
    Int top_n_splicing = 25
    Int top_n_protein = 15
  }

  call PrepareWindowGenotypes {
    input:
      window_id = window_id,
      trans_window_associations = trans_window_associations,
      genome_dosage = genome_dosage,
      genome_dosage_tbi = genome_dosage_tbi
  }

  call PrepareWindowPhenotypes {
    input:
      window_id = window_id,
      trans_window_associations = trans_window_associations,
      expression_phenotypes = expression_phenotypes,
      splicing_phenotypes = splicing_phenotypes,
      protein_phenotypes = protein_phenotypes,
      target_phenotypes = target_phenotypes,
      expression_phenotypes_tbi = expression_phenotypes_tbi,
      expression_phenotype_lookup = expression_phenotype_lookup,
      splicing_phenotypes_tbi = splicing_phenotypes_tbi,
      splicing_phenotype_lookup = splicing_phenotype_lookup,
      protein_phenotypes_tbi = protein_phenotypes_tbi,
      protein_phenotype_lookup = protein_phenotype_lookup,
      top_n_expression = top_n_expression,
      top_n_splicing = top_n_splicing,
      top_n_protein = top_n_protein
  }

  output {
    File window_dosage = PrepareWindowGenotypes.window_dosage
    File window_manifest = PrepareWindowGenotypes.window_manifest
    File window_phenotypes = PrepareWindowPhenotypes.window_phenotypes
    File phenotype_data = PrepareWindowPhenotypes.phenotype_data
    File window_qc = PrepareWindowPhenotypes.window_qc
  }
}

task PrepareWindowGenotypes {
  input {
    File trans_window_associations
    String window_id
    File genome_dosage
    File genome_dosage_tbi
  }

  command <<<
    set -euo pipefail

    output_prefix='~{window_id}'
    if [[ ! "$output_prefix" =~ ^[A-Za-z0-9._-]+$ ]]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: window_id contains unsafe filename characters." >&2
      exit 1
    fi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting genotype preparation for ~{window_id}."
    mkdir -p output

    dosage_name="window_input.dose.tsv.gz"
    ln -sf "~{genome_dosage}" "${dosage_name}"
    ln -sf "~{genome_dosage_tbi}" "${dosage_name}.tbi"
    test -s "${dosage_name}"
    test -s "${dosage_name}.tbi"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Local dosage: ${dosage_name}."
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Local index: ${dosage_name}.tbi."

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Resolving the locus coordinates."
    window_row="$(awk -F '\t' -v requested_id='~{window_id}' '
      NR == 1 {
        for (i = 1; i <= NF; i++) column[$i] = i
        next
      }
      $(column["window_id"]) == requested_id {
        row = $(column["chrom"]) "\t" $(column["start"]) "\t" $(column["end"])
        if (!seen[row]++) {
          selected = row
          matches++
        }
      }
      END {
        if (matches != 1) exit 1
        print selected
      }
    ' <(gzip -cd ~{trans_window_associations}))"

    IFS=$'\t' read -r window_chrom window_start window_end <<< "${window_row}"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Interval: ${window_chrom}:$((window_start + 1))-${window_end}."
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Reading the dosage header."
    tabix -H "${dosage_name}" > dosage_header.raw.tsv
    if [[ -s dosage_header.raw.tsv ]]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using the header returned by tabix."
    else
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] tabix returned no header. Reading the first compressed line."
      set +o pipefail
      gzip -cd -- "${dosage_name}" \
        | awk 'NR == 1 {print; exit}' \
        > dosage_header.raw.tsv
      set -o pipefail
    fi
    test -s dosage_header.raw.tsv
    header_line="$(tail -n 1 dosage_header.raw.tsv)"
    header_line="${header_line#\#}"
    printf '%s\n' "${header_line}" > "output/${output_prefix}.window_dosage.tsv"
    if ! awk -F '\t' '
      NR == 1 {
        valid = ($1 == "CHROM" && $2 == "POS" && $3 == "REF" && $4 == "ALT")
        exit !valid
      }
    ' "output/${output_prefix}.window_dosage.tsv"; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Dosage header must start with CHROM, POS, REF, and ALT." >&2
      exit 1
    fi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Dosage header validation complete."

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracting all locus genotype rows."
    tabix "${dosage_name}" \
      "${window_chrom}:$((window_start + 1))-${window_end}" \
      >> "output/${output_prefix}.window_dosage.tsv"
    test "$(wc -l < "output/${output_prefix}.window_dosage.tsv")" -gt 1
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Extracted variants: $(( $(wc -l < "output/${output_prefix}.window_dosage.tsv") - 1 ))."

    {
      printf 'window_id\tchrom\tstart\tend\tdosage_file\n'
      printf '%s\t%s\t%s\t%s\t%s\n' \
        '~{window_id}' \
        "${window_chrom}" \
        "${window_start}" \
        "${window_end}" \
        "${output_prefix}.window_dosage.tsv"
    } > "output/${output_prefix}.window_manifest.tsv"
    test -s "output/${output_prefix}.window_manifest.tsv"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Genotype preparation complete."
  >>>

  output {
    File window_dosage = "output/" + window_id + ".window_dosage.tsv"
    File window_manifest = "output/" + window_id + ".window_manifest.tsv"
  }

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/mvsusier-prepare-window-genotypes:latest"
    cpu: 2
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
  }
}

task PrepareWindowPhenotypes {
  input {
    String window_id
    File trans_window_associations
    File expression_phenotypes
    File splicing_phenotypes
    File protein_phenotypes
    File target_phenotypes
    File? expression_phenotypes_tbi
    File? expression_phenotype_lookup
    File? splicing_phenotypes_tbi
    File? splicing_phenotype_lookup
    File? protein_phenotypes_tbi
    File? protein_phenotype_lookup
    Int top_n_expression
    Int top_n_splicing
    Int top_n_protein
  }

  command <<<
    set -euo pipefail

    output_prefix='~{window_id}'
    if [[ ! "$output_prefix" =~ ^[A-Za-z0-9._-]+$ ]]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: window_id contains unsafe filename characters." >&2
      exit 1
    fi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting joint phenotype preparation for ~{window_id}."
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Top-N values: expression=~{top_n_expression}, splicing=~{top_n_splicing}, protein=~{top_n_protein}."

    expression_tbi='~{default="" expression_phenotypes_tbi}'
    expression_lookup='~{default="" expression_phenotype_lookup}'
    splicing_tbi='~{default="" splicing_phenotypes_tbi}'
    splicing_lookup='~{default="" splicing_phenotype_lookup}'
    protein_tbi='~{default="" protein_phenotypes_tbi}'
    protein_lookup='~{default="" protein_phenotype_lookup}'

    access_mode() {
      if [[ -n "$1" && -n "$2" ]]; then
        printf 'tabix'
      elif [[ -z "$1" && -z "$2" ]]; then
        printf 'full_scan'
      else
        printf 'incomplete'
      fi
    }
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Phenotype access: expression=$(access_mode "$expression_tbi" "$expression_lookup"), splicing=$(access_mode "$splicing_tbi" "$splicing_lookup"), protein=$(access_mode "$protein_tbi" "$protein_lookup")."

    optional_args=()
    if [[ -n "$expression_tbi" ]]; then
      optional_args+=(--expression-phenotypes-tbi "$expression_tbi")
    fi
    if [[ -n "$expression_lookup" ]]; then
      optional_args+=(--expression-phenotype-lookup "$expression_lookup")
    fi
    if [[ -n "$splicing_tbi" ]]; then
      optional_args+=(--splicing-phenotypes-tbi "$splicing_tbi")
    fi
    if [[ -n "$splicing_lookup" ]]; then
      optional_args+=(--splicing-phenotype-lookup "$splicing_lookup")
    fi
    if [[ -n "$protein_tbi" ]]; then
      optional_args+=(--protein-phenotypes-tbi "$protein_tbi")
    fi
    if [[ -n "$protein_lookup" ]]; then
      optional_args+=(--protein-phenotype-lookup "$protein_lookup")
    fi

    Rscript /opt/mvsusie/scripts/prepare_trans_window.R \
      --window-id ~{window_id} \
      --trans-associations ~{trans_window_associations} \
      --expression-phenotypes ~{expression_phenotypes} \
      --splicing-phenotypes ~{splicing_phenotypes} \
      --protein-phenotypes ~{protein_phenotypes} \
      --target-phenotypes ~{target_phenotypes} \
      --top-n-expression ~{top_n_expression} \
      --top-n-splicing ~{top_n_splicing} \
      --top-n-protein ~{top_n_protein} \
      "${optional_args[@]}" \
      --output-dir output
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Validating joint phenotype outputs."
    test -s "output/${output_prefix}.window_phenotypes.tsv"
    test -s "output/${output_prefix}.window_phenotypes.bed.gz"
    test -s "output/${output_prefix}.window_qc.tsv"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Joint phenotype preparation complete."
  >>>

  output {
    File window_phenotypes = "output/" + window_id + ".window_phenotypes.tsv"
    File phenotype_data = "output/" + window_id + ".window_phenotypes.bed.gz"
    File window_qc = "output/" + window_id + ".window_qc.tsv"
  }

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/mvsusier-prepare-window-phenotypes:latest"
    cpu: 2
    memory: "16 GiB"
    disks: "local-disk 500 SSD"
  }
}
