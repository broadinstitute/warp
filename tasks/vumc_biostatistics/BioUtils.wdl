version 1.0

task GetGeneLocus {
  input {
    String gene_symbol

    String docker = "shengqh/r4:20241117"
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
    File output_bed = "~{target_file}"
  }
}
