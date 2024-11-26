version 1.0

task RegenieStep1FitModel {
  input {
    File input_pgen
    File input_pvar
    File input_psam

    File phenoFile
    String phenoColList
    Boolean is_binary_traits

    File covarFile
    String covarColList

    # Regenie options
    # option "--loocv" is not in the recommendation of Regenie (https://rgcgithub.github.io/regenie/recommendations/)
    # However, using --loocv would accelerate the process a lot. We still suggest to use it.
    #############################################################
    # Options in effect:
    # --step 1 \
    # --qt \
    # --pgen /cromwell_root/fc-9b4e856a-12ef-40a9-aca0-8d5f9c2ab9c1/submissions/470ae9a0-7258-49c1-a1f5-c2281d1a2855/VUMCRegenie/45f954cc-a4da-4dfd-8601-599fe18d48e1/call-PgenQCFilter/cacheCopy/demo_bmi_953.qc \
    # --phenoFile /cromwell_root/fc-0a566538-7eb6-45c9-94a8-c4002fcb63ce/demo/gwas/demo_bmi_953_phenotype.final.txt \
    # --phenoColList bmi \
    # --covarFile /cromwell_root/fc-0a566538-7eb6-45c9-94a8-c4002fcb63ce/demo/gwas/demo_bmi_953_phenotype.final.txt \
    # --covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    # --bsize 1000 \
    # --lowmem \
    # --threads 8 \
    # --out demo_bmi_953.step1 \
    # --force-step1
    # Chromosome 1
    # block [1] : 1000 snps (5ms)
    # -residualizing and scaling genotypes...done (9ms)
    # -calc working matrices...done (56ms)
    # -calc level 0 ridge...done (543ms)
    #############################################################
    # Options in effect:
    #   --step 1 \
    #   --qt \
    #   --pgen /cromwell_root/fc-9b4e856a-12ef-40a9-aca0-8d5f9c2ab9c1/submissions/6094d580-cc50-4811-9fb7-cedd399970c2/VUMCRegenie/71318aa2-90b9-49c1-880c-7ca522b5f86b/call-PgenQCFilter/cacheCopy/demo_bmi_953.qc \
    #   --phenoFile /cromwell_root/fc-0a566538-7eb6-45c9-94a8-c4002fcb63ce/demo/gwas/demo_bmi_953_phenotype.final.txt \
    #   --phenoColList bmi \
    #   --covarFile /cromwell_root/fc-0a566538-7eb6-45c9-94a8-c4002fcb63ce/demo/gwas/demo_bmi_953_phenotype.final.txt \
    #   --covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    #   --loocv \
    #   --bsize 1000 \
    #   --lowmem \
    #   --threads 8 \
    #   --out demo_bmi_953.step1 \
    #   --force-step1
    # Chromosome 1
    #  block [1] : 1000 snps  (5ms) 
    #    -residualizing and scaling genotypes...done (9ms) 
    #    -calc working matrices...done (135ms) 
    #    -calc level 0 ridge...done (18ms) 

    String step1_option = "--loocv --bsize 1000 --lowmem"

    String output_prefix

    Int memory_gb = 40
    Int cpu = 8

    String docker = "skoyamamd/regenie:3.4.2"
  }

  Int disk_size = ceil(size([input_pgen, input_pvar, input_psam], "GB")) + 10

  String call_type = if(is_binary_traits) then "--bt" else "--qt"

  command <<<

pgen='~{input_pgen}'
pgen_prefix=${pgen%.*}

regenie --step 1 \
  ~{call_type} \
  --pgen ${pgen_prefix} \
  -p ~{phenoFile} \
  --phenoColList ~{phenoColList} \
  -c ~{covarFile} \
  --covarColList ~{covarColList} \
  ~{step1_option} \
  --threads ~{cpu} \
  --out ~{output_prefix} \
  --force-step1

>>>

  runtime {
    docker: docker
    preemptible: 1
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  output {
    File pred_list_file = "${output_prefix}_pred.list" 
    Array[File] pred_loco_files = glob("*.loco")
  }
}

task RegenieStep2AssociationTest {
  input {
    File input_pgen
    File input_pvar
    File input_psam

    File pred_list_file
    Array[File] pred_loco_files

    File phenoFile
    String phenoColList
    Boolean is_binary_traits

    File covarFile
    String covarColList

    String step2_option = "--firth --approx --pThresh 0.01 --bsize 400"

    Int? chromosome

    String output_prefix

    Int memory_gb = 100
    Int cpu = 8

    String docker = "skoyamamd/regenie:3.4.2"
  }

  Int disk_size = ceil(size([input_pgen, input_pvar, input_psam], "GB")) + 20

  String call_type = if(is_binary_traits) then "--bt" else "--qt"

  command <<<

set -euo pipefail
        
for file in ~{sep=' ' pred_loco_files}; do \
  mv $file .; \
done

while IFS=' ' read -r name path; do
  filename=$(basename "$path")
  echo "$name $filename"
done < ~{pred_list_file} > pred.list

pgen='~{input_pgen}'
pgen_prefix=${pgen%.*}
 
regenie --step 2 \
  ~{call_type} ~{"--chr " + chromosome} \
  --pgen ${pgen_prefix} \
  -p ~{phenoFile} \
  --phenoColList ~{phenoColList} \
  -c ~{covarFile} \
  --covarColList ~{covarColList} \
  ~{step2_option} \
  --threads ~{cpu} \
  --pred pred.list \
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
    Array[File] regenie_files = glob("~{output_prefix}*.regenie")
  }
}
