version 1.0

task GetGeneLocus {
  input {
    String gene_symbol
    Int shift_bases = 2000

    String docker = "shengqh/report:20241120"
    Int preemptible = 1

    String host = "https://www.ensembl.org"
    String dataset = "hsapiens_gene_ensembl"
    String symbolKey = "hgnc_symbol"
    Int addChr = 1
  }

  String target_file = gene_symbol + ".bed"

  command <<<

cat <<EOF > script.r

require(biomaRt)
require(stringr)

host="~{host}"
dataset="~{dataset}"
symbolKey="~{symbolKey}"
genes="~{gene_symbol}"
addChr=~{addChr}
shift_bases=~{shift_bases}

ensembl <- useMart("ensembl", host=host, dataset=dataset)

geneLocus<-getBM(attributes=c("chromosome_name", "start_position", "end_position", symbolKey, "strand", "ensembl_gene_id"),
                 filters=symbolKey, 
                 values=genes, 
                 mart=ensembl, 
                 uniqueRows=TRUE,
                 useCache=FALSE)

geneLocus<-geneLocus[nchar(geneLocus\$chromosome_name) < 6,]

geneLocus\$score<-1000

geneLocus<-geneLocus[,c("chromosome_name", "start_position", "end_position", "score", symbolKey, "strand", "ensembl_gene_id")]
geneLocus<-geneLocus[order(geneLocus\$chromosome_name, geneLocus\$start_position),]

geneLocus\$strand[geneLocus\$strand == 1]<-"+"
geneLocus\$strand[geneLocus\$strand == -1]<-"-"

if(addChr & (!any(grepl("chr", geneLocus\$chromosome_name)))){
  geneLocus\$chromosome_name = paste0("chr", geneLocus\$chromosome_name)
}

if(shift_bases > 0){
  geneLocus\$start_position = geneLocus\$start_position - shift_bases
  geneLocus\$end_position = geneLocus\$end_position + shift_bases
}

bedFile<-"~{target_file}"
write.table(geneLocus, file=bedFile, row.names=F, col.names = F, sep="\t", quote=F)

EOF

R -f script.r

>>>

  runtime {
    cpu: 1
    docker: "~{docker}"
    preemptible: preemptible
    disks: "local-disk 10 HDD"
    memory: "4 GiB"
  }
  output {
    File gene_bed = "~{target_file}"
  }
}

task PgenQCFilter {
  input {
    File input_pgen
    File input_pvar
    File input_psam

    String output_prefix

    String qc_option

    Int memory_gb = 20
    Int cpu = 8
    Float disk_size_factor = 1.5

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([input_pgen, input_pvar, input_psam], "GB")  * disk_size_factor) + 20

  command <<<

plink2 \
  --pgen ~{input_pgen} \
  --pvar ~{input_pvar} \
  --psam ~{input_psam} \
  ~{qc_option} \
  --make-pgen \
  --out ~{output_prefix}

>>>

  runtime {
    docker: docker
    preemptible: 1
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File output_pgen = "~{output_prefix}.pgen"
    File output_pvar = "~{output_prefix}.pvar"
    File output_psam = "~{output_prefix}.psam"
  }
}
