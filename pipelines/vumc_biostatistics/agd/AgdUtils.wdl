version 1.0

task ReplaceICAIdWithGrid {
  input {
    File input_psam
    File id_map_file
    String output_psam
  }

  command <<<

python3 <<CODE

import io

with open("~{id_map_file}", "rt") as fin:
  id_map = {}
  for line in fin:
    parts = line.strip().split('\t')
    id_map[parts[0]] = parts[1]

with open("~{input_psam}", "rt") as fin:
  with open("~{output_psam}", "wt") as fout:
    for line in fin:
      parts = line.strip().split('\t')
      if parts[1] in id_map:
        grid = id_map[parts[1]]
        parts[1] = grid

      newline = '\t'.join(parts)
      fout.write(f"{newline}\n")

CODE

>>>

  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    File output_psam = "~{output_psam}"
  }
}

task CreateCohortPsam {
  input {
    File input_psam
    File grid_file
    String output_psam
  }

  command <<<

python3 <<CODE

import os

grids = set(line.strip() for line in open("~{grid_file}", "rt"))
with open("~{output_psam}", "wt") as fout:
  with open("~{input_psam}", "rt") as fin:
    for line in fin:
      if line.startswith("#"):
        fout.write(line)
        continue
      if line.split('\t')[1] in grids:
        fout.write(line)
CODE

echo "Number of samples to keep:"
wc -l "~{output_psam}"

>>>

  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    preemptible: 1
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
  }
  output {
    File output_psam = "~{output_psam}"
  }
}


task HailMatrixExtractRegions {
  input {
    File input_hail_mt_path_file
    File input_bed
    Float expect_output_vcf_bgz_size_gb
    String target_prefix

    String billing_project_id

    String docker = "shengqh/hail_gcp:20240211"
    Int memory_gb = 10
    Int preemptible = 1
    Int cpu = 4
    Int boot_disk_gb = 10  
  }

  Int disk_size = ceil(expect_output_vcf_bgz_size_gb) + 20
  Int total_memory_gb = memory_gb + 2

  String target_file = target_prefix + ".vcf.bgz"

  command <<<

#https://discuss.hail.is/t/i-get-a-negativearraysizeexception-when-loading-a-plink-file/899
export PYSPARK_SUBMIT_ARGS="--driver-java-options '-XX:hashCode=0' --conf 'spark.executor.extraJavaOptions=-XX:hashCode=0' pyspark-shell"

python3 <<CODE

import hail as hl
import pandas as pd

def parse_gcs_url(gcs_url):
    if not gcs_url.startswith('gs://'):
        raise ValueError("URL must start with 'gs://'")

    # Remove the 'gs://' prefix
    gcs_url = gcs_url[5:]

    # Split the remaining URL into bucket name and object key
    parts = gcs_url.split('/', 1)
    if len(parts) != 2:
        raise ValueError("Invalid GCS URL format")

    bucket_name = parts[0]

    return bucket_name

regions = pd.read_csv("~{input_bed}", sep='\t', header=None)
regions=regions.iloc[:, :3]
regions.columns = ["chr", "start", "end"]
regions['chr']=regions['chr'].astype(str)
regions['locus']=regions.chr + ":" + (regions.start + 1).astype(str) + "-" + (regions.end + 1).astype(str)
regions.head()

new_tbl = pd.read_csv("~{input_hail_mt_path_file}", sep='\t')
new_tbl.head()

hail_url = new_tbl['hail'][0]
if hail_url.startswith('gs://'):
  bucket_name = parse_gcs_url(new_tbl['hail'][0])
  print(f"hail_bucket_name={bucket_name}")

  hl.init(spark_conf={"spark.driver.memory": "~{memory_gb}g",
                      'spark.hadoop.fs.gs.requester.pays.mode': 'CUSTOM',
                      'spark.hadoop.fs.gs.requester.pays.buckets': bucket_name,
                      'spark.hadoop.fs.gs.requester.pays.project.id': "~{billing_project_id}"}, idempotent=True)
else:
  hl.init(spark_conf={"spark.driver.memory": "~{memory_gb}g"}, idempotent=True)

hl.default_reference("GRCh38")

hail_col="hail"

all_tbl=None
for ind in new_tbl.index:
    chr=new_tbl['chromosome'][ind]
    chr_regions=regions[regions.chr==chr]
    if chr_regions.shape[0] == 0:
        print(f"{chr}: no snps")
    else:
        print(f"{chr}: {chr_regions.shape[0]} snps")
        print(chr_regions)
        hail_url=new_tbl[hail_col][ind]
        mt = hl.read_matrix_table(hail_url)

        mt_filter = hl.filter_intervals(
            mt,
            [hl.parse_locus_interval(x,) for x in chr_regions.locus])
        
        print(f"  {chr} found {mt_filter.count_rows()} snps from hailmatrix")

        if mt_filter.count_rows() > 0:
            if all_tbl == None:
                all_tbl = mt_filter
            else:
                all_tbl = all_tbl.union_rows(mt_filter)

print(f"SNPs={all_tbl.count_rows()}, Samples={all_tbl.count_cols()}")

print(f"Writing to ~{target_file}")

hl.export_vcf(all_tbl, "~{target_file}")

CODE

>>>

  runtime {
    cpu: cpu
    docker: "~{docker}"
    preemptible: preemptible
    disks: "local-disk ~{disk_size} HDD"
    memory: "~{total_memory_gb} GiB"
    bootDiskSizeGb: boot_disk_gb
  }
  output {
    File output_vcf = "~{target_file}"
  }
}

task Annovar {
  input {
    File input_vcf
    String annovar_db = "/opt/annovar/humandb"
    String annovar_option = ""

    String target_prefix

    Int memory_gb = 20
    Int cpu = 8

    String docker = "shengqh/annovar:20241117"
  }

  Int disk_size = ceil(size([input_vcf], "GB")  * 2) + 10

  command <<<

zcat ~{input_vcf} | cut -f1-9 > ~{target_prefix}.avinput.vcf

convert2annovar.pl -format vcf4old ~{target_prefix}.avinput.vcf | cut -f1-7 | awk '{gsub(",\\*", "", $0); print}'> ~{target_prefix}.avinput

table_annovar.pl ~{target_prefix}.avinput ~{annovar_db} -buildver hg38 -protocol refGene,avsnp150 -operation g,f --remove  --outfile ~{target_prefix}.annovar

>>>

  runtime {
    docker: docker
    preemptible: 1
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File output_annovar_file = "~{target_prefix}.annovar.hg38_multianno.txt"
  }
}

