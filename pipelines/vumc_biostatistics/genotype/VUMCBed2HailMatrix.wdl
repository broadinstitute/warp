version 1.0

workflow VUMCBed2HailMatrix {
  input {
    File source_bed
    File source_bim
    File source_fam

    String reference_genome = "GRCh38"

    String target_prefix

    String? project_id
    String target_gcp_folder
  }

  call Bed2HailMatrix {
    input:
      source_bed = source_bed,
      source_bim = source_bim,
      source_fam = source_fam,

      reference_genome = reference_genome,
      target_prefix = target_prefix,

      project_id = project_id,
      target_gcp_folder = target_gcp_folder
  }

  output {
    String hail_gcs_path = Bed2HailMatrix.hail_gcs_path
  }
}

task Bed2HailMatrix {
  input {
    File source_bed
    File source_bim
    File source_fam

    String reference_genome
    String target_prefix

    String? project_id
    String target_gcp_folder

    String docker = "hailgenetics/hail:0.2.127-py3.11"
    Int memory_gb = 20
  }

  Int disk_size = ceil(size([source_bed, source_bim, source_fam], "GB")  * 3) + 20
  Int total_memory_gb = memory_gb + 2

  String gcs_output_dir = sub(target_gcp_folder, "/+$", "")
  String gcs_output_path = gcs_output_dir + "/" + target_prefix

  command <<<

#https://discuss.hail.is/t/i-get-a-negativearraysizeexception-when-loading-a-plink-file/899

export PYSPARK_SUBMIT_ARGS="--driver-java-options '-XX:hashCode=0' --conf 'spark.executor.extraJavaOptions=-XX:hashCode=0' pyspark-shell"

python3 <<CODE

import hail as hl

hl.init(spark_conf={"spark.driver.memory": "~{memory_gb}g"})

#contig_recoding is hard coded for human only
dsplink = hl.import_plink(bed="~{source_bed}",
                          bim="~{source_bim}",
                          fam="~{source_fam}",
                          reference_genome="~{reference_genome}",
                          contig_recoding={
                            '1': 'chr1',
                            '2': 'chr2',
                            '3': 'chr3',
                            '4': 'chr4',
                            '5': 'chr5',
                            '6': 'chr6',
                            '7': 'chr7',
                            '8': 'chr8',
                            '9': 'chr9',
                            '10': 'chr10',
                            '11': 'chr11',
                            '12': 'chr12',
                            '13': 'chr13',
                            '14': 'chr14',
                            '15': 'chr15',
                            '16': 'chr16',
                            '17': 'chr17',
                            '18': 'chr18',
                            '19': 'chr19',
                            '20': 'chr20',
                            '21': 'chr21',
                            '22': 'chr22',
                            'X': 'chrX',
                            'Y': 'chrY',
                            'MT': 'chrM'})

dsplink.write("~{target_prefix}", overwrite=True)

CODE

gsutil ~{"-u " + project_id} -m rsync -Cr ~{target_prefix} ~{gcs_output_path}

#tar czf ~{target_prefix}.tar.gz ~{target_prefix} 

>>>

  runtime {
    docker: "~{docker}"
    preemptible: 1
    disks: "local-disk ~{disk_size} HDD"
    memory: "~{total_memory_gb} GiB"
  }
  output {
    String hail_gcs_path = "~{gcs_output_path}"
  }
}
