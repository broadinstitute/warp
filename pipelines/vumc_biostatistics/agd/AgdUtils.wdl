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
    String? target_gcp_folder

    String docker = "shengqh/hail_gcp:20240211"
    Int memory_gb = 10
    Int preemptible = 1
    Int cpu = 4
    Int boot_disk_gb = 10  
  }

  Int disk_size = ceil(expect_output_vcf_bgz_size_gb * 2) + 20
  Int total_memory_gb = memory_gb + 2

  String target_file = if defined(target_gcp_folder) then sub(select_first([target_gcp_folder]), "/+$", "") + "/" + target_prefix + ".vcf.bgz" else target_prefix + ".vcf.bgz"

  command <<<

#https://discuss.hail.is/t/i-get-a-negativearraysizeexception-when-loading-a-plink-file/899
export PYSPARK_SUBMIT_ARGS="--driver-java-options '-XX:hashCode=0' --conf 'spark.executor.extraJavaOptions=-XX:hashCode=0' pyspark-shell"

mkdir tmp

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

  hl.init(tmp_dir='tmp',
          spark_conf={"spark.driver.memory": "~{memory_gb}g",
                      "spark.local.dir": "tmp",
                      'spark.hadoop.fs.gs.requester.pays.mode': 'CUSTOM',
                      'spark.hadoop.fs.gs.requester.pays.buckets': bucket_name,
                      'spark.hadoop.fs.gs.requester.pays.project.id': "~{billing_project_id}"}, idempotent=True)
else:
  hl.init(tmp_dir='tmp',
          spark_conf={"spark.driver.memory": "~{memory_gb}g",
                      "spark.local.dir": "tmp"}, idempotent=True)

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
            #keep GT only
            mt_filter = mt_filter.select_entries(mt_filter.GT)

            # Remove INFO annotations by setting them to empty
            mt_filter = mt_filter.annotate_rows(info=hl.struct())

            if all_tbl == None:
                all_tbl = mt_filter
            else:
                all_tbl = all_tbl.union_rows(mt_filter)

print(f"SNPs={all_tbl.count_rows()}, Samples={all_tbl.count_cols()}")

print(f"Writing to ~{target_file}")

all_tbl = all_tbl.naive_coalesce(2)

hl.export_vcf(all_tbl, "~{target_file}")

CODE

# #keep GT only
# bcftools annotate -x INFO,^FORMAT/GT ~{target_file} -Ov -o tmp.vcf.bgz
# mv tmp.vcf.bgz ~{target_file}

rm -rf tmp

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
    Int cpu = 1

    String docker = "shengqh/annovar:20241117"
  }

  Int disk_size = ceil(size([input_vcf], "GB")  * 2) + 10

  command <<<

zcat ~{input_vcf} | cut -f1-9 > ~{target_prefix}.avinput.vcf

convert2annovar.pl -format vcf4old ~{target_prefix}.avinput.vcf | cut -f1-7 | awk '{gsub(",\\*", "", $0); print}'> ~{target_prefix}.avinput

table_annovar.pl ~{target_prefix}.avinput ~{annovar_db} -buildver hg38 -protocol refGene -operation g --remove  --outfile ~{target_prefix}.annovar

>>>

  runtime {
    docker: docker
    preemptible: 1
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File annovar_file = "~{target_prefix}.annovar.hg38_multianno.txt"
  }
}

task PrepareGeneGenotype {
  input {
    String gene_symbol
    File agd_primary_grid_file
    File annovar_file
    File vcf_file

    String docker = "shengqh/report:20241118"
    
    Int preemptible = 1
    Int additional_disk_size = 10
    Int memory_gb = 20
  }

  Int disk_size = ceil(size([agd_primary_grid_file, annovar_file, vcf_file], "GB") * 5) + additional_disk_size

  command <<<

cat <<EOF > script.r

library(data.table)
library(dplyr)

gene='~{gene_symbol}'
agd_primary_grid_file='~{agd_primary_grid_file}'
annovar_file='~{annovar_file}'
vcf_file='~{vcf_file}'

cat("reading", agd_primary_grid_file, "..\n")
agd_df=fread(agd_primary_grid_file, data.table=FALSE)

cat("reading", annovar_file, "...\n")
annovar=fread(annovar_file)
cat("there are", nrow(annovar), "SNVs in annovar ...\n")

print(head(annovar))

cat("filtering by gene", gene, "...\n")
annovar = annovar |>
  dplyr::filter(Gene.refGene==gene)
cat("there are", nrow(annovar), "SNVs in annovar from", gene, "...\n")

cat("filtering snv ... \n")
lof_snv = rbind(annovar |> dplyr::filter(Func.refGene %in% c('splicing')),
            annovar |> dplyr::filter(Func.refGene %in% c('exonic')) |> dplyr::filter(ExonicFunc.refGene %in% c('stopgain', 'startloss'))
)
vuc_snv = rbind(annovar |> dplyr::filter(Func.refGene %in% c('splicing')),
            annovar |> dplyr::filter(Func.refGene %in% c('exonic')) |> dplyr::filter(ExonicFunc.refGene %in% c('stopgain', 'startloss', 'nonsynonymous SNV'))
)

cat("there are", nrow(vuc_snv), "valid vuc SNVs, including", nrow(lof_snv), "valid lof SNVs.\n")

# Use a pipe to decompress with zcat and read the first 4000 lines
con <- pipe(paste("zcat", vcf_file, "| head -n 4000"), "rt")
first_lines <- readLines(con)

# Close the connection
close(con)

chrom_index=grep("^#CHROM", first_lines)
cat("data starts from line", chrom_index, "...\n")

cat("reading", vcf_file, "...\n")
vcf = fread(cmd=paste0("zcat ", vcf_file), skip=chrom_index-1, data.table=FALSE)
cat("there are total", nrow(vcf), "SNVs...\n")

to_genotype_file<-function(vcf, snv, genotype_file){
  cat("preparing", genotype_file, "...\n")
  snv_vcf=vcf |> dplyr::filter(POS %in% snv\$Start)
  snv_vcf_data = snv_vcf[,10:ncol(snv_vcf)]

  cat("  converting snv to genotype ... \n")
  snv_vcf_gt = data.frame(lapply(snv_vcf_data, function(x) { gsub(':.*', '', x)}), check.names=FALSE)
  snv_vcf_gt = data.frame(lapply(snv_vcf_data, function(x) { gsub('[|]', '/', x)}), check.names=FALSE)
  print(head(snv_vcf_gt[,1:5]))

  has_snv=apply(snv_vcf_gt, 2, function(x) { any(x %in% c('1/1', '0/1', '1/0', '0/2', '2/0'))})

  df=data.frame(GRID=colnames(snv_vcf_gt), Genotype=ifelse(has_snv, "1", "0")) |> 
    dplyr::filter(GRID %in% agd_df\$PRIMARY_GRID) |>
    dplyr::arrange(GRID)
  print(table(df\$Genotype))
    
  cat("  saving to", genotype_file, "...\n")
  write.csv(df, genotype_file, quote=FALSE, row.names=FALSE)
}

to_genotype_file(vcf, lof_snv, paste0(gene, ".lof.genotype.csv"))
to_genotype_file(vcf, vuc_snv, paste0(gene, ".vuc.genotype.csv"))

cat("done\n")

EOF

R -f script.r

>>>

  runtime {
    cpu: 1
    docker: "~{docker}"
    preemptible: preemptible
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File lof_genotype_file = "~{gene_symbol}.lof.genotype.csv"
    File vuc_genotype_file = "~{gene_symbol}.vuc.genotype.csv"
  }
}

task PreparePhenotype {
  input {
    String phename
    Float phecode

    File agd_primary_grid_file
    File phecode_data_file
    File phecode_map_file
    Int min_occurance = 2

    String docker = "shengqh/report:20241118"
    
    Int preemptible = 1

    Int memory_gb = 20
  }

  Int disk_size = ceil(size([agd_primary_grid_file, phecode_data_file, phecode_map_file], "GB")) + 10

  command <<<

wget https://raw.githubusercontent.com/shengqh/ngsperl/refs/heads/master/lib/CQS/reportFunctions.R

wget https://raw.githubusercontent.com/shengqh/ngsperl/refs/heads/master/lib/BioVU/prepare_phenotype_data.rmd

mv prepare_phenotype_data.rmd ~{phename}.phenotype.rmd

echo -e "~{phename}\tphename" > input_options.txt
echo -e "~{phecode}\tphecode" >> input_options.txt
echo -e "~{agd_primary_grid_file}\tagd_file" >> input_options.txt
echo -e "~{phecode_data_file}\tphecode_data_file" >> input_options.txt
echo -e "~{phecode_map_file}\tphecode_map_file" >> input_options.txt
echo -e "~{min_occurance}\tmin_occurance" >> input_options.txt

R -e "library(knitr);rmarkdown::render(input='~{phename}.phenotype.rmd');"   

>>>

  runtime {
    cpu: 1
    docker: "~{docker}"
    preemptible: preemptible
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File phenotype_file = "~{phename}.phenotype.csv"
    File phenotype_report = "~{phename}.phenotype.html"
  }
}

task LinearAssociation {
  input {
    String phename
    Float phecode
    String genotype_name

    File phenotype_file

    File agd_primary_grid_file
    File demographics_file
    File genotype_file
    File pca_file
    File phecode_map_file
    File ancestry_file

    String docker = "shengqh/report:20241118"
    
    Int preemptible = 1

    Int memory_gb = 20
  }

  Int disk_size = ceil(size([phenotype_file, agd_primary_grid_file, demographics_file, genotype_file, pca_file, phecode_map_file, ancestry_file], "GB")) + 10

  command <<<

wget https://raw.githubusercontent.com/shengqh/ngsperl/refs/heads/master/lib/CQS/reportFunctions.R

wget https://raw.githubusercontent.com/shengqh/ngsperl/refs/heads/master/lib/BioVU/linear_association.rmd

mv linear_association.rmd ~{phecode}.~{genotype_name}.glm.rmd

echo -e "~{phename}\tphename" >> input_options.txt
echo -e "~{phecode}\tphecode" > input_options.txt
echo -e "~{phenotype_file}\tphefile" >> input_options.txt
echo -e "~{genotype_name}\tgenotype_name" >> input_options.txt
echo -e "~{genotype_file}\tgenotype_file" >> input_options.txt
echo -e "~{agd_primary_grid_file}\tagd_file" >> input_options.txt
echo -e "~{demographics_file}\tdemographics_file" >> input_options.txt
echo -e "~{pca_file}\tpca_file" >> input_options.txt
echo -e "~{phecode_map_file}\tphecode_map_file" >> input_options.txt
echo -e "~{ancestry_file}\tancestry_file" >> input_options.txt

R -e "library(knitr);rmarkdown::render(input='~{phename}.~{genotype_name}.glm.rmd');"   

>>>

  runtime {
    cpu: 1
    docker: "~{docker}"
    preemptible: preemptible
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File linear_association_file = "~{phename}.~{genotype_name}.glm.csv"
    File linear_association_report = "~{phename}.~{genotype_name}.glm.html"
  }
}
