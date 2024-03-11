version 1.0

import "./Utils.wdl" as Utils

workflow VUMCVcf2GDS {
  input {
    File source_vcf
    File source_vcf_index

    File? replace_samples

    String reference_genome = "GRCh38"

    String target_prefix

    String target_suffix = ".gds"

    String docker = "shengqh/gds:20240205"

    Int memory_gb = 20

    String? project_id
    String? target_gcp_folder
  }

  call Vcf2GDS {
    input:
      source_vcf = source_vcf,
      source_vcf_index = source_vcf_index,
      replace_samples = replace_samples,
      reference_genome = reference_genome,
      target_prefix = target_prefix,
      target_suffix = target_suffix,
      docker = docker,
      memory_gb = memory_gb
  }

  if(defined(target_gcp_folder)){
    call Utils.MoveOrCopyTwoFiles as CopyFile {
      input:
        source_file1 = Vcf2GDS.gds_file,
        source_file2 = Vcf2GDS.gds_samples,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File gds_file = select_first([CopyFile.output_file1, Vcf2GDS.gds_file])
    File gds_samples = select_first([CopyFile.output_file2, Vcf2GDS.gds_samples])
  }
}

task Vcf2GDS {
  input {
    File source_vcf
    File source_vcf_index

    File? replace_samples

    String reference_genome
    String target_prefix
    String target_suffix
    String docker
    Int memory_gb
    Int cpu = 8
  }

  Int disk_size = ceil(size([source_vcf, source_vcf_index], "GB")  * 1.5) + 20
  Int total_memory_gb = memory_gb + 2
  String target_file = "~{target_prefix}~{target_suffix}"

  command <<<

cat <<EOT >> convertVcftoGDS.r

library(SeqArray)

seqVCF2GDS("~{source_vcf}", "~{target_file}", reference="~{reference_genome}", parallel=~{cpu})

if ("~{replace_samples}" != "") {
  gds <- seqOpen("~{target_file}")
  old_sample_ids = seqGetData(gds, "sample.id")
  seqClose(gds)

  writeLines(old_sample_ids, "~{target_file}.old_samples.txt")

  idmap=read.table("~{replace_samples}", sep="\t", header=FALSE)
  idmap\$V1=as.character(idmap\$V1)
  idmap\$V2=as.character(idmap\$V2)
  if(!all(old_sample_ids %in% idmap\$V1)){
    stop("Not all samples in the GDS file are in the replace_samples file")
  }

  new_sample_ids = idmap\$V2[match(old_sample_ids, idmap\$V1)]

  f <- openfn.gds("~{target_file}", readonly=FALSE)
  add.gdsn(f, "sample.id", new_sample_ids, replace=TRUE, compress="LZMA_RA", closezip=TRUE)
  closefn.gds(f)
}

gds <- seqOpen("~{target_file}")
new_sample_ids = seqGetData(gds, "sample.id")
seqClose(gds)

writeLines(new_sample_ids, "~{target_file}.samples.txt")

EOT

R --vanilla -f convertVcftoGDS.r

>>>

  runtime {
    docker: "~{docker}"
    cpu: cpu
    preemptible: 1
    disks: "local-disk ~{disk_size} HDD"
    memory: "~{total_memory_gb} GiB"
  }
  output {
    File gds_file = "~{target_file}"
    File gds_samples = "~{target_file}.samples.txt"
  }
}
