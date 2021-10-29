version 1.0

# Tasks that are expected to be used only on Broad Infrastructure

# This task is to convert a string into another string which is 'safe' as a filename.
# That is, it does not contain any number of odd characters.
# Note that we are not testing for the occurrence of a single quote (') in the string.
task MakeSafeFilename {
  input {
    String name
  }

  String safe_name = name

  command <<<
    echo '~{name}' | tr ' "#$%&*/:;<=>?@[]^{}|~\\()' '_' > safe_name.txt
  >>>

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 10 HDD"
    memory: "1 GiB"
    preemptible: 3
  }
  output {
    String output_safe_name = read_string('safe_name.txt')
  }
}

task DownloadGenotypes {
  input {
    String sample_alias
    String sample_lsid
    String output_vcf_base_name
    File params_file
    Boolean compress = true

    File haplotype_database_file

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String environment
    File vault_token_path

    Int? max_retries
    Int? preemptible_tries
  }

  meta {
    volatile: true
  }

  String fp_retrieved_file = "fp_retrieved.txt"

  String output_vcf = output_vcf_base_name + if compress then ".vcf.gz" else".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"

  command <<<

    export VAULT_ADDR=https://clotho.broadinstitute.org:8200
    export VAULT_TOKEN=$(cat ~{vault_token_path})
    if [ ~{environment} == prod ]; then
      export MERCURY_AUTH_KEY=secret/dsde/gotc/prod/wdl/secrets
      export MERCURY_FP_STORE_URI=https://portals.broadinstitute.org/portal/mercury-ws/fingerprint
    else
      export MERCURY_AUTH_KEY=secret/dsde/gotc/dev/wdl/secrets
      export MERCURY_FP_STORE_URI=https://portals.broadinstitute.org/portal-test/mercury-ws/fingerprint
    fi

    grep 'REGULATORY_DESIGNATION=RESEARCH_ONLY' ~{params_file}
    if [ $? -eq 0 ]; then
      ADDITIONAL_PLATFORM_STRING="--EXPECTED_GENOTYPING_PLATFORMS GENERAL_ARRAY"
    fi

    exit_code=0

    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
    DownloadGenotypes \
      --SAMPLE_ALIAS "~{sample_alias}" \
      --SAMPLE_LSID "~{sample_lsid}" \
      --OUTPUT "~{output_vcf}" \
      --CREATE_INDEX true \
      --REFERENCE_SEQUENCE ~{ref_fasta} \
      --HAPLOTYPE_MAP ~{haplotype_database_file} \
      --EXPECTED_GENOTYPING_PLATFORMS FLUIDIGM \
      $ADDITIONAL_PLATFORM_STRING \
      --IGNORE_SPECIFIC_GENOTYPES_PLATFORM GENERAL_ARRAY \
      --IGNORE_SPECIFIC_GENOTYPES_LSID ~{sample_lsid} \
      --MERCURY_FP_STORE_URI $MERCURY_FP_STORE_URI \
      --CREDENTIALS_VAULT_PATH $MERCURY_AUTH_KEY \
      --ERR_NO_GENOTYPES_AVAILABLE 7
    exit_code=$?

    if [ $exit_code -eq 7 ]; then
      # Exit code from DownloadGenotypes if no fingerprints were found.
      # Treat this as a normal condition, but set a variable to indicate no fingerprints available.
      # Create empty file so that it exists.
      exit_code=0
      echo "Found no fingerprints for ~{sample_lsid}!"
      echo "false" > ~{fp_retrieved_file}
      touch ~{output_vcf}
      touch ~{output_vcf_index}
    else
      echo "true" > ~{fp_retrieved_file}
    fi

    exit $exit_code
  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.1.0-1631191359"
    memory: "3.5 GiB"
    maxRetries: select_first([max_retries, 2])
    preemptible: select_first([preemptible_tries, 3])
  }

  output {
    Boolean fingerprint_retrieved = read_boolean(fp_retrieved_file)
    File reference_fingerprint_vcf = output_vcf
    File reference_fingerprint_vcf_index = output_vcf_index
  }
}


task UploadFingerprintToMercury {
  input {
    File fingerprint_json_file
    File gtc_file

    String environment
    File vault_token_path

    Int? max_retries
    Int? preemptible_tries
  }

  command <<<
    set -eo pipefail

    export VAULT_ADDR=https://clotho.broadinstitute.org:8200
    export VAULT_TOKEN=$(cat ~{vault_token_path})
    if [ ~{environment} == prod ]; then
      export MERCURY_AUTH_KEY=secret/dsde/gotc/prod/wdl/secrets
      export MERCURY_FP_STORE_URI=https://portals.broadinstitute.org/portal/mercury-ws/fingerprint
    else
      export MERCURY_AUTH_KEY=secret/dsde/gotc/dev/wdl/secrets
      export MERCURY_FP_STORE_URI=https://portals.broadinstitute.org/portal-test/mercury-ws/fingerprint
    fi

    du -k ~{gtc_file} | cut -f 1 > size.txt

    # TODO -Fix UploadFingerprintToMercury so I don't need to pass a file size

    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/gitc/picard-private.jar \
      UploadFingerprintToMercury \
      --INPUT "~{fingerprint_json_file}" \
      --GTC_FILE_SIZE size.txt \
      --MERCURY_FP_STORE_URI $MERCURY_FP_STORE_URI \
      --CREDENTIALS_VAULT_PATH $MERCURY_AUTH_KEY \
  >>>

  runtime {
    docker: "us.gcr.io/broad-arrays-prod/arrays-picard-private:4.1.0-1631191359"
    memory: "3.5 GiB"
    maxRetries: select_first([max_retries, 2])
    preemptible: select_first([preemptible_tries, 3])
  }
}
