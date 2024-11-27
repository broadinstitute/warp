version 1.0

import "../../../tasks/vumc_biostatistics/GcpUtils.wdl" as GcpUtils

workflow VUMCAGDPhecode {
  input {
    String bigquery_project_id
    String bigquery_dataset_id
    String id_mapping_file
    File phemap_file
    String output_prefix

    String? billing_gcp_project_id
    String? target_gcp_folder    
  }

  call query_phecode {
    input:
      bigquery_project_id = bigquery_project_id,
      bigquery_dataset_id = bigquery_dataset_id,
      id_mapping_file = id_mapping_file,
      phemap_file = phemap_file,
      output_prefix = output_prefix
  }

  if(defined(target_gcp_folder)){
    String gcs_output_dir = select_first([target_gcp_folder])

    call GcpUtils.MoveOrCopyFourFiles as CopyFile {
      input:
        source_file1 = query_phecode.demographics_csv,
        source_file2 = query_phecode.phecode12_long_csv,
        source_file3 = query_phecode.phecode12_wide_csv,
        source_file4 = query_phecode.phecode12_wide_binarized_csv,
        is_move_file = false,
        project_id = billing_gcp_project_id,
        target_gcp_folder = gcs_output_dir
    }
  }

  output {
    File demographics_csv = select_first([CopyFile.output_file1, query_phecode.demographics_csv])
    File phecode12_long_csv = select_first([CopyFile.output_file2, query_phecode.phecode12_long_csv])
    File phecode12_wide_csv = select_first([CopyFile.output_file3, query_phecode.phecode12_wide_csv])
    File phecode12_wide_binarized_csv = select_first([CopyFile.output_file4, query_phecode.phecode12_wide_binarized_csv])
  }
}

task query_phecode {
  input {
    String bigquery_project_id
    String bigquery_dataset_id
    
    String id_mapping_file
    File phemap_file

    String output_prefix

    String docker = "shengqh/hail_gcp:20241127"
  }

  command <<<

cat <<EOF > query.py

import pandas as pd
import pandas_gbq

from google.cloud import bigquery

# Configure the external data source and query job.
external_config = bigquery.ExternalConfig("CSV")
external_config.source_uris = [
    "~{id_mapping_file}"
]
external_config.schema = [
    bigquery.SchemaField("GRID", "STRING"), 
    bigquery.SchemaField("sample_id","NUMERIC")
]
assert external_config.csv_options is not None
external_config.csv_options.skip_leading_rows = 1

table_id = "terra_temp"
job_config = bigquery.QueryJobConfig(table_definitions={table_id: external_config})

client = bigquery.Client()

sql = 'SELECT * FROM {}'.format(table_id)
print(sql)

results = client.query(sql, job_config=job_config).result().to_dataframe()
print(results.head())

#demographics for the grids in that dataset
sd_url="~{bigquery_project_id}.~{bigquery_dataset_id}"
query_sql=f"""select *
            from
            {sd_url}.person as p
            inner join terra_temp as tt
                on p.person_source_value=tt.GRID;"""

res = client.query(query_sql, job_config=job_config).result().to_dataframe() 
res.to_csv("~{output_prefix}.demographics.csv",index=False)

phemap = pd.read_csv("~{phemap_file}",sep='\t',dtype=str)
phemap.head()

#get icd codes for the mapping file set
icd_sql=f"""select person_source_value, concept_code as icd, vocabulary_id
    from 
    {sd_url}.all_events as ae
    inner join {sd_url}.concept as c
        on c.concept_id=ae.source_concept_id
    inner join {sd_url}.person as p
        on p.person_id=ae.person_id
    inner join terra_temp as tt
        on p.person_source_value=tt.GRID
    where c.vocabulary_id like 'ICD%CM';"""

icd_codes = client.query(icd_sql, job_config=job_config).result().to_dataframe() 
icd_codes.head()

phemap['vocabulary_id']='ICD'+phemap['flag'].astype(str)+'CM'
merged_old = icd_codes.merge(phemap, on=['ICD','vocabulary_id'], how='inner')
merged_old.head()

wide = pd.pivot_table(merged_old[['GRID','phecode']], index='GRID', columns = 'phecode', aggfunc='size',fill_value=0)
transformed = wide.applymap(transform_counts)

merged_old[['GRID','phecode']].to_csv("~{output_prefix}.Phecode12_long.csv",index=False)
wide.to_csv("~{output_prefix}.Phecode12_wide.csv")
transformed.to_csv("~{output_prefix}.Phecode12_wide_binarized.csv")

EOF

python3 query.py

  >>>

  runtime {
    docker: docker
    memory: 10 + " GiB"
    disks: "local-disk " + 10 + " HDD"
    cpu: 1
    preemptible: 1
  }

  output {
    File demographics_csv = "~{output_prefix}.demographics.csv"
    File phecode12_long_csv = "~{output_prefix}.Phecode12_long.csv"
    File phecode12_wide_csv = "~{output_prefix}.Phecode12_wide.csv"
    File phecode12_wide_binarized_csv = "~{output_prefix}.Phecode12_wide_binarized.csv"
  }
}
