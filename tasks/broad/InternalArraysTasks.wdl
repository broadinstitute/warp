version 1.0

task CreateExtendedIlluminaManifest {
  input {
    File input_csv
    String output_base_name
    File cluster_file

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbSNP_vcf_file
    File dbSNP_vcf_index_file

    File supported_ref_fasta
    File supported_ref_fasta_index
    File supported_ref_dict

    File chain_file

    Int disk_size
    Int preemptible_tries
  }

  command <<<
    java -Xms13g -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
            CreateExtendedIlluminaManifest \
            --INPUT ~{input_csv} \
            --OUTPUT_BASE_FILE ~{output_base_name} \
            --CLUSTER_FILE ~{cluster_file} \
            --DBSNP_FILE ~{dbSNP_vcf_file} \
            --TARGET_BUILD 37 \
            --TARGET_REFERENCE_FILE ~{ref_fasta} \
            --SUPPORTED_BUILD 36 \
            --SUPPORTED_REFERENCE_FILE ~{supported_ref_fasta} \
            --SUPPORTED_CHAIN_FILE ~{chain_file}
  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.0.10-1602016912"
    disks: "local-disk " + disk_size + " HDD"
    memory: "14 GiB"
    preemptible: preemptible_tries
  }

  output {
    File output_csv = "~{output_base_name}.extended.csv"
    File output_bad_assays_file = "~{output_base_name}.bad_assays.csv"
    File report_file = "~{output_base_name}.report.txt"
  }
}

task GenerateEmptyVariantCallingMetricsFile {
  input {
    String chip_well_barcode
    String sample_alias
    String chip_type
    String reported_gender
    String autocall_version
    String output_metrics_basename
    String cluster_filename
    Int analysis_version_number
    Int preemptible_tries
  }

  command <<<
    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
            GenerateEmptyVariantCallingMetrics \
            --CHIP_WELL_BARCODE ~{chip_well_barcode} \
            --SAMPLE_ALIAS "~{sample_alias}" \
            --CHIP_TYPE ~{chip_type} \
            --REPORTED_GENDER "~{reported_gender}" \
            --CLUSTER_FILE_NAME ~{cluster_filename} \
            --AUTOCALL_VERSION ~{autocall_version} \
            --ANALYSIS_VERSION ~{analysis_version_number} \
            --OUTPUT ~{output_metrics_basename}
  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.0.10-1602016912"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
      File detail_metrics = "~{output_metrics_basename}.arrays_variant_calling_detail_metrics"
    }
}

task BlacklistBarcode {
  input {
    File upload_metrics_output
    String chip_well_barcode
    Int analysis_version
    Int preemptible_tries
    Array[String] authentication
    String service_account_filename
    String reason
    String notes
  }

  command <<<
    set -eo pipefail

    AUTH=~{write_lines(authentication)} && chmod +x $AUTH && $AUTH
    export GOOGLE_APPLICATION_CREDENTIALS=/cromwell_root/~{service_account_filename}

    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
                  ArraysManualBlacklistUpdate \
                  --CHIP_WELL_BARCODE ~{chip_well_barcode} \
                  --ANALYSIS_VERSION ~{analysis_version} \
                  --REASON ~{reason} \
                  --DB_USERNAME_FILE cloudsql.db_user.txt \
                  --DB_PASSWORD_FILE cloudsql.db_password.txt \
                  --DB_JDBC_FILE cloudsql.db_jdbc.txt \
                  --NOTES "~{notes}"
  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.0.10-1602016912"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }
}

task VcfToMercuryFingerprintJson {
  input {
    File input_vcf_file
    File input_vcf_index_file
    File variant_calling_detail_metrics_file
    String sample_lsid
    String output_json_filename
    Int disk_size
    Int preemptible_tries
  }

  command <<<
    set -eo pipefail

    # Need to determine the disposition from the metrics
    # Remove all the comments and blank lines from the file
    # Find the column number of AUTOCALL_PF and retrieve the value in the second line of the AUTOCALL_PF column
    # AUTOCALL_PF set to empty if file has more than 2 lines (should only have column headers and one data line)
    AUTOCALL_PF=$(sed '/#/d' ~{variant_calling_detail_metrics_file} |
      sed '/^\s*$/d' |
      awk -v col=AUTOCALL_PF -F '\t' \
      'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}}} NR==2{print $c} NR>2{exit 1}')

    if [[ "$AUTOCALL_PF" == "Y" ]]
    then
       DISPOSITION=P
    elif [[ "$AUTOCALL_PF" == "N" ]]
    then
       DISPOSITION=F
    else
       echo "AUTOCALL_PF should only be Y or N and there should only be one line of data"
       exit 1;
    fi

    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
      VcfToMercuryFingerprintJson \
      --VCF ~{input_vcf_file} \
      --LSID ~{sample_lsid} \
      --PLATFORM GENERAL_ARRAY \
      --DISPOSITION $DISPOSITION \
      --SNP_LIST_NAME HG19HaplotypeDbSnps \
      --OUTPUT ~{output_json_filename}
  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.0.10-1602016912"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
    File output_json_file = output_json_filename
  }
}


task CreateBafRegressMetricsFile {
  input {
    File input_file
    String output_metrics_basefilename

    Int disk_size
    Int preemptible_tries
  }

  command {
    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
      CreateBafRegressMetricsFile \
      --INPUT ~{input_file} \
      --OUTPUT ~{output_metrics_basefilename}
  }
  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.0.10-1602016912"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
    File output_metrics_file = "~{output_metrics_basefilename}.bafregress_metrics"
  }
}

task UploadArraysMetrics {
  input {
    File arrays_variant_calling_detail_metrics
    File? arrays_variant_calling_summary_metrics
    File? arrays_control_code_summary_metrics
    File? fingerprinting_detail_metrics
    File? fingerprinting_summary_metrics
    File? genotype_concordance_summary_metrics
    File? genotype_concordance_detail_metrics
    File? genotype_concordance_contingency_metrics
    File? verify_id_metrics
    File? bafregress_metrics

    Array[String] authentication
    String service_account_filename

    Int disk_size
    Int preemptible_tries
  }

  command <<<
    set -eo pipefail

    AUTH=~{write_lines(authentication)} && chmod +x $AUTH && $AUTH
    export GOOGLE_APPLICATION_CREDENTIALS=/cromwell_root/~{service_account_filename}

    rm -rf metrics_upload_dir &&
    mkdir metrics_upload_dir &&

    # check that files are passed in before copying them -- [ -z FILE ] evaluates to true if FILE not there
    ! [ -z ~{genotype_concordance_summary_metrics} ] &&
    cp ~{genotype_concordance_summary_metrics} metrics_upload_dir
    ! [ -z ~{genotype_concordance_detail_metrics} ] &&
    cp ~{genotype_concordance_detail_metrics} metrics_upload_dir
    ! [ -z ~{genotype_concordance_contingency_metrics} ] &&
    cp ~{genotype_concordance_contingency_metrics} metrics_upload_dir
    ! [ -z ~{verify_id_metrics} ] &&
    cp ~{verify_id_metrics} metrics_upload_dir
    ! [ -z ~{bafregress_metrics} ] &&
    cp ~{bafregress_metrics} metrics_upload_dir

    ! [ -z ~{fingerprinting_detail_metrics} ] &&
    cp ~{fingerprinting_detail_metrics} metrics_upload_dir
    ! [ -z ~{fingerprinting_summary_metrics} ] &&
    cp ~{fingerprinting_summary_metrics} metrics_upload_dir

    cp ~{arrays_variant_calling_detail_metrics} metrics_upload_dir
    ! [ -z ~{arrays_variant_calling_summary_metrics} ] &&
    cp ~{arrays_variant_calling_summary_metrics} metrics_upload_dir

    ! [ -z ~{arrays_control_code_summary_metrics} ] &&
    cp ~{arrays_control_code_summary_metrics} metrics_upload_dir
    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
      UploadArraysMetrics \
      --ANALYSIS_DIRECTORY metrics_upload_dir \
      --DB_USERNAME_FILE cloudsql.db_user.txt \
      --DB_PASSWORD_FILE cloudsql.db_password.txt \
      --DB_JDBC_FILE cloudsql.db_jdbc.txt &&
    touch empty_file_for_dependency
  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.0.10-1602016912"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
      File upload_metrics_empty_file = "empty_file_for_dependency"
    }
}

task UpdateChipWellBarcodeIndex {
  input {
    File params_file
    Array[String] authentication
    String service_account_filename
    Int disk_size
    Int preemptible_tries
  }

  command <<<
    set -eo pipefail

    AUTH=~{write_lines(authentication)} && chmod +x $AUTH && $AUTH
    export GOOGLE_APPLICATION_CREDENTIALS=/cromwell_root/~{service_account_filename}
    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
      UpdateChipWellBarcodeIndex \
      --PARAMS_FILE ~{params_file} \
      --DB_USERNAME_FILE cloudsql.db_user.txt \
      --DB_PASSWORD_FILE cloudsql.db_password.txt \
      --DB_JDBC_FILE cloudsql.db_jdbc.txt

  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.0.10-1602016912"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }
}
