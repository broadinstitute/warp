version 1.0

workflow VUMCBigQueryDemo {
  input {
    String project_id='sd-vumc-tanagra-test'
    String dataset_id='terra_sd_20240831'
    String table_id='person'
  }

  call query_table_from_bigquery{
    input:
      project_id = project_id,
      dataset_id = dataset_id,
      table_id = table_id
  }

  output {
    File query_result_json = query_table_from_bigquery.query_result_json
  }
}

task query_table_from_bigquery {
  input {
    String project_id
    String dataset_id
    String table_id

    String docker = "google/cloud-sdk:latest"
  }

  command {
    # The `bq` command will automatically authenticate using the service account 
    # associated with your Terra workspace, so no need for explicit credentials setup.
    bq query --nouse_legacy_sql "SELECT * FROM \`~{project_id}.~{dataset_id}.~{table_id}\` LIMIT 10" > result.json
  }

  runtime {
    docker: docker
    memory: 10 + " GiB"
    disks: "local-disk " + 10 + " HDD"
    cpu: 1
    preemptible: 1
  }

  output {
    File query_result_json = "result.json"
  }
}
