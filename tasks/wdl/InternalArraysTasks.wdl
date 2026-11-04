version 1.0

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
    java -Xms2000m -Xmx3000m -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
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
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.1.3-1652895718"
    memory: "3500 MiB"
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
    Int analysis_version_number
    Int preemptible_tries
    File vault_token_path
    Array[String] authentication
    String service_account_filename
    String reason
    String notes
  }

  meta {
    volatile: true
  }

  command <<<
    set -eo pipefail

    export VAULT_TOKEN=$(cat ~{vault_token_path})
    AUTH=~{write_lines(authentication)} && chmod +x $AUTH && $AUTH
    export GOOGLE_APPLICATION_CREDENTIALS=/cromwell_root/~{service_account_filename}

    java -Xms2000m -Xmx3000m -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
                  ArraysManualBlacklistUpdate \
                  --CHIP_WELL_BARCODE ~{chip_well_barcode} \
                  --ANALYSIS_VERSION ~{analysis_version_number} \
                  --REASON ~{reason} \
                  --DB_USERNAME_FILE cloudsql.db_user.txt \
                  --DB_PASSWORD_FILE cloudsql.db_password.txt \
                  --DB_JDBC_FILE cloudsql.db_jdbc.txt \
                  --NOTES "~{notes}"
  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.1.3-1652895718"
    memory: "3500 MiB"
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

    java -Xms2000m -Xmx3000m -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
      VcfToMercuryFingerprintJson \
      --VCF ~{input_vcf_file} \
      --LSID ~{sample_lsid} \
      --PLATFORM GENERAL_ARRAY \
      --DISPOSITION $DISPOSITION \
      --SNP_LIST_NAME HG19HaplotypeDbSnps \
      --OUTPUT ~{output_json_filename}
  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.1.3-1652895718"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3500 MiB"
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
    java -Xms2000m -Xmx3000m -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
      CreateBafRegressMetricsFile \
      --INPUT ~{input_file} \
      --OUTPUT ~{output_metrics_basefilename}
  }
  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.1.3-1652895718"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3500 MiB"
    preemptible: preemptible_tries
  }

  output {
    File output_metrics_file = "~{output_metrics_basefilename}.bafregress_metrics"
  }
}

task UploadArraysMetrics {
  input {
    File arrays_variant_calling_detail_metrics
    File arrays_variant_calling_summary_metrics
    File arrays_control_code_summary_metrics
    File? fingerprinting_detail_metrics
    File? fingerprinting_summary_metrics
    File? genotype_concordance_summary_metrics
    File? genotype_concordance_detail_metrics
    File? genotype_concordance_contingency_metrics
    File? verify_id_metrics
    File? bafregress_metrics

    File vault_token_path
    Array[String] authentication
    String service_account_filename

    Int disk_size
    Int preemptible_tries
  }

  meta {
    volatile: true
  }

  command <<<
    set -eo pipefail

    export VAULT_TOKEN=$(cat ~{vault_token_path})
    AUTH=~{write_lines(authentication)} && chmod +x $AUTH && $AUTH
    export GOOGLE_APPLICATION_CREDENTIALS=/cromwell_root/~{service_account_filename}

    rm -rf metrics_upload_dir &&
    mkdir metrics_upload_dir &&

    cp ~{arrays_control_code_summary_metrics} metrics_upload_dir
    cp ~{arrays_variant_calling_detail_metrics} metrics_upload_dir
    cp ~{arrays_variant_calling_summary_metrics} metrics_upload_dir

    # check that optional files exist before copying them -- [ -z FILE ] evaluates to true if FILE not there
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

    java -Xms2000m -Xmx3000m -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
      UploadArraysMetrics \
      --ANALYSIS_DIRECTORY metrics_upload_dir \
      --DB_USERNAME_FILE cloudsql.db_user.txt \
      --DB_PASSWORD_FILE cloudsql.db_password.txt \
      --DB_JDBC_FILE cloudsql.db_jdbc.txt &&
    touch empty_file_for_dependency
  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.1.3-1652895718"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3500 MiB"
    preemptible: preemptible_tries
  }

  output {
      File upload_metrics_empty_file = "empty_file_for_dependency"
    }
}

task UploadEmptyArraysMetrics {
  input {
    File arrays_variant_calling_detail_metrics

    File vault_token_path
    Array[String] authentication
    String service_account_filename

    Int disk_size
    Int preemptible_tries
  }

  meta {
    volatile: true
  }

  command <<<
    set -eo pipefail

    export VAULT_TOKEN=$(cat ~{vault_token_path})
    AUTH=~{write_lines(authentication)} && chmod +x $AUTH && $AUTH
    export GOOGLE_APPLICATION_CREDENTIALS=/cromwell_root/~{service_account_filename}

    rm -rf metrics_upload_dir &&
    mkdir metrics_upload_dir &&

    cp ~{arrays_variant_calling_detail_metrics} metrics_upload_dir

    java -Xms2000m -Xmx3000m -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
    UploadArraysMetrics \
    --ANALYSIS_DIRECTORY metrics_upload_dir \
    --DB_USERNAME_FILE cloudsql.db_user.txt \
    --DB_PASSWORD_FILE cloudsql.db_password.txt \
    --DB_JDBC_FILE cloudsql.db_jdbc.txt &&
    touch empty_file_for_dependency
  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.1.3-1652895718"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3500 MiB"
    preemptible: preemptible_tries
  }

  output {
    File upload_metrics_empty_file = "empty_file_for_dependency"
  }
}

task CreateChipWellBarcodeParamsFile {
  input {
    String chip_type_name
    String chip_well_barcode
    String? collaborator_participant_id
    String? lab_batch
    String? participant_id
    String? product_family
    String? product_name
    String? product_order_id
    String? product_part_number
    String product_type
    String regulatory_designation
    String research_project_id
    String sample_alias
    String gender
    String? sample_id
    String sample_lsid
    Int preemptible_tries
  }

  String params_filename = "params.txt"

  command <<<
    set -eo pipefail

    echo "CHIP_TYPE_NAME=~{chip_type_name}" > ~{params_filename}
    echo "CHIP_WELL_BARCODE=~{chip_well_barcode}" >> ~{params_filename}
    echo "INDIVIDUAL_ALIAS=~{collaborator_participant_id}" >> ~{params_filename}
    echo "LAB_BATCH=~{lab_batch}" >> ~{params_filename}
    echo "PARTICIPANT_ID=~{participant_id}" >> ~{params_filename}
    echo "PRODUCT_FAMILY=~{product_family}" >> ~{params_filename}
    echo "PRODUCT_NAME=~{product_name}" >> ~{params_filename}
    echo "PRODUCT_ORDER_ID=~{product_order_id}" >> ~{params_filename}
    echo "PRODUCT_PART_NUMBER=~{product_part_number}" >> ~{params_filename}
    echo "PRODUCT_TYPE=~{product_type}" >> ~{params_filename}
    echo "REGULATORY_DESIGNATION=~{regulatory_designation}" >> ~{params_filename}
    echo "RESEARCH_PROJECT_ID=~{research_project_id}" >> ~{params_filename}
    echo "SAMPLE_ALIAS=~{sample_alias}" >> ~{params_filename}
    echo "SAMPLE_GENDER=~{gender}" >> ~{params_filename}
    echo "SAMPLE_ID=~{sample_id}" >> ~{params_filename}
    echo "SAMPLE_LSID=~{sample_lsid}" >> ~{params_filename}

  >>>

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: preemptible_tries
  }
  output {
    File params_file = params_filename
  }
}



task UpdateChipWellBarcodeIndex {
  input {
    File params_file
    File vault_token_path
    Array[String] authentication
    String service_account_filename
    Int disk_size
    Int preemptible_tries
  }

  meta {
    volatile: true
  }

  command <<<
    set -eo pipefail

    export VAULT_TOKEN=$(cat ~{vault_token_path})
    AUTH=~{write_lines(authentication)} && chmod +x $AUTH && $AUTH
    export GOOGLE_APPLICATION_CREDENTIALS=/cromwell_root/~{service_account_filename}
    java -Xms2000m -Xmx3000m -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
      UpdateChipWellBarcodeIndex \
      --PARAMS_FILE ~{params_file} \
      --DB_USERNAME_FILE cloudsql.db_user.txt \
      --DB_PASSWORD_FILE cloudsql.db_password.txt \
      --DB_JDBC_FILE cloudsql.db_jdbc.txt

  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.1.3-1652895718"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3500 MiB"
    preemptible: preemptible_tries
  }
}

task GetNextArraysQcAnalysisVersionNumber {
  input {
    String chip_well_barcode
    Int preemptible_tries
    File vault_token_path
    Array[String] authentication
    String service_account_filename
  }

  meta {
    volatile: true
  }

  command <<<
    set -eo pipefail

    export VAULT_TOKEN=$(cat ~{vault_token_path})
    AUTH=~{write_lines(authentication)} && chmod +x $AUTH && $AUTH
    export GOOGLE_APPLICATION_CREDENTIALS=/cromwell_root/~{service_account_filename}

    java -Xms2000m -Xmx3000m -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
      GetNextArraysQcAnalysisVersionNumber \
        --CHIP_WELL_BARCODE ~{chip_well_barcode} \
        --DB_USERNAME_FILE cloudsql.db_user.txt \
        --DB_PASSWORD_FILE cloudsql.db_password.txt \
        --DB_JDBC_FILE cloudsql.db_jdbc.txt
  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.1.3-1652895718"
    memory: "3500 MiB"
    preemptible: preemptible_tries
  }
  output {
    Int analysis_version_number = read_int(stdout())
  }
}

task ResolveExtendedIlluminaManifestFile {
  input {
    File extended_manifest_map_file
    String bpm_filename
    String egt_filename
    String arrays_chip_metadata_path
    Int preemptible_tries
  }

  command <<<
    cat ~{extended_manifest_map_file} | grep -v '^#' | grep -E '^~{bpm_filename}\s+~{egt_filename}' | cut -f 3 > output_file.txt
    if [[ ! -s output_file.txt ]]
    then
      echo "ERROR: Unable to find entry in ~{extended_manifest_map_file} for ~{bpm_filename} / ~{egt_filename}" 1>&2
      exit 1
    elif [[ $(cat output_file.txt | wc -l) -ne 1 ]]
    then
      echo "ERROR: Found more than one entry in ~{extended_manifest_map_file} for ~{bpm_filename} / ~{egt_filename}" 1>&2
      exit 1
    fi
  >>>

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: preemptible_tries
  }
  output {
    String extended_illumina_manifest_file = arrays_chip_metadata_path + read_string("output_file.txt")
  }
}

task ResolveMinorAlleleFrequencyFile {
  input {
    File minor_allele_frequency_map_file
    String bpm_filename
    String arrays_chip_metadata_path
    Int preemptible_tries
  }

  command <<<
    cat ~{minor_allele_frequency_map_file} | grep -v '^#' | grep ~{bpm_filename} | cut -f 2 > output_file.txt
    if [[ ! -s output_file.txt ]]
    then
      echo "Unable to find entry in ~{minor_allele_frequency_map_file} for ~{bpm_filename}"
      echo false > found.txt
      exit 0
    elif [[ $(cat output_file.txt | wc -l) -ne 1 ]]
    then
      echo "ERROR: Found more than one entry in ~{minor_allele_frequency_map_file} for ~{bpm_filename}" 1>&2
      echo false > found.txt
      exit 1
    fi
    echo true > found.txt
  >>>

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: preemptible_tries
  }
  output {
    Boolean found = read_boolean("found.txt")
    String minor_allele_frequency_file = arrays_chip_metadata_path + read_string("output_file.txt")
  }
}

task FormatArraysOutputs {
    input {
        String  chip_well_barcode_output
        Int     analysis_version_number_output
        String? baf_regress_metrics_file
        String? gtc_file

        String? output_vcf
        String? output_vcf_index

        String? arrays_variant_calling_detail_metrics_file
        String? arrays_variant_calling_summary_metrics_file
        String? arrays_variant_calling_control_metrics_file

        String? fingerprint_detail_metrics_file
        String? fingerprint_summary_metrics_file

        String? genotype_concordance_summary_metrics_file
        String? genotype_concordance_detail_metrics_file
        String? genotype_concordance_contingency_metrics_file

        String  lab_batch
    }

    command <<<
        echo -e "chip_well_barcode_output\tanalysis_version_number_output\tlab_batch\tbaf_regress_metrics_file\tgtc_file\t\
        output_vcf\toutput_vcf_index\t\
        arrays_variant_calling_detail_metrics_file\tarrays_variant_calling_summary_metrics_file\tarrays_variant_calling_control_metrics_file\t\
        fingerprint_detail_metrics_file\tfingerprint_summary_metrics_file\t\
        genotype_concordance_summary_metrics_file\tgenotype_concordance_detail_metrics_file\tgenotype_concordance_contingency_metrics_file" \
        > ingestDataset_arrays_outputs.tsv

        echo -e "~{chip_well_barcode_output}\t~{analysis_version_number_output}\t~{lab_batch}\t~{baf_regress_metrics_file}\t~{gtc_file}\t\
        ~{output_vcf}\t~{output_vcf_index}\t\
        ~{arrays_variant_calling_detail_metrics_file}\t~{arrays_variant_calling_summary_metrics_file}\t~{arrays_variant_calling_control_metrics_file}\t\
        ~{fingerprint_detail_metrics_file}\t~{fingerprint_summary_metrics_file}\t\
        ~{genotype_concordance_summary_metrics_file}\t~{genotype_concordance_detail_metrics_file}\t~{genotype_concordance_contingency_metrics_file}" \
        >> ingestDataset_arrays_outputs.tsv

        python3 << CODE
        import pandas as pd

        tsv_df = pd.read_csv("ingestDataset_arrays_outputs.tsv", sep="\t")
        tsv_df = tsv_df.dropna(axis=1, how="all")  # drop columns if no value (optional outputs etc)

        outputs = tsv_df.to_json("ingestDataset_arrays_outputs.json", orient="records")  # write json file

        CODE
    >>>

    runtime {
        docker: "gcr.io/emerge-production/emerge_wdls:v.1.0"
    }

    output {
        File ingest_outputs_tsv = "ingestDataset_arrays_outputs.tsv"
        File ingest_outputs_json = "ingestDataset_arrays_outputs.json"
    }
}