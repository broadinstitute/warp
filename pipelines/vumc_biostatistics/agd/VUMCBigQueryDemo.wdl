version 1.0

workflow VUMCBigQueryDemo {
  input {
    String project_id='sd-vumc-tanagra-test'
    String dataset_id='terra_sd_20240831'
    String table_id='person'
  }

  call query_table_by_bq {
    input:
      project_id = project_id,
      dataset_id = dataset_id,
      table_id = table_id
  }

  call query_table_by_pandas {
    input:
      project_id = project_id,
      dataset_id = dataset_id,
      table_id = table_id
  }

  output {
    File query_result_json = query_table_by_bq.query_result_json
    File query_result_csv = query_table_by_pandas.query_result_csv
  }
}

task query_table_by_bq {
  input {
    String project_id
    String dataset_id
    String table_id

    String docker = "shengqh/hail_gcp:20241120"
  }

  command {
    # The `bq` command will automatically authenticate using the service account 
    # associated with your Terra workspace, so no need for explicit credentials setup.
    bq query --nouse_legacy_sql "SELECT * FROM \`~{project_id}.~{dataset_id}.~{table_id}\` LIMIT 1" > result.json
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

task query_table_by_pandas {
  input {
    String project_id
    String dataset_id
    String table_id

    String docker = "shengqh/hail_gcp:20241120"
  }

  command <<<

cat <<EOF > query_table.py

import pandas as pd
from google.cloud import bigquery

client = bigquery.Client()
sql=f"SELECT * FROM `~{project_id}.~{dataset_id}.~{table_id}` LIMIT 1"
patient_data=client.query(sql).result().to_dataframe() 
patient_data.to_csv('result.csv', index=False)

EOF

python query_table.py

  >>>

  runtime {
    docker: docker
    memory: 10 + " GiB"
    disks: "local-disk " + 10 + " HDD"
    cpu: 1
    preemptible: 1
  }

  output {
    File query_result_csv = "result.csv"
  }
}
