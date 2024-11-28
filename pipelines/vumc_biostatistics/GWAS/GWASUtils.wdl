version 1.0

task Regenie4MemoryEstimation {
  input {
    Int num_sample
    Int num_variant
    Int num_phenotype
    Int num_chromosome
    Int num_covariate
    Int num_ridge_l0 = 5
    Int block_size = 1000

    String docker = "shengqh/report:20241120"
  }

  command <<<
    cat <<EOF > script.r

num_sample=~{num_sample}
num_variant=~{num_variant}
num_phenotype=~{num_phenotype}
num_chromosome=~{num_chromosome}
num_covariate=~{num_covariate}
num_ridge_l0=~{num_ridge_l0}
block_size=~{block_size}

# https://github.com/rgcgithub/regenie/blob/5bb86f9a207422e4affafb9eb17693c87ec30c0a/src/Regenie.cpp
# // Step 1
# // 4P + max( B + PRT, PRT) + #chrs [P:#traits;R=#ridge l0;T=#predictions from l0]
# int t_eff = ( params->write_l0_pred ? 1 : params->total_n_block );

total_n_block = ceiling(1.0 * num_variant / block_size)
t_eff = total_n_block

# int p_eff = ( params->write_l0_pred ? 1 : params->n_pheno );

p_eff = num_phenotype

# int b_eff = params->total_n_block;

b_eff = total_n_block

# total_ram = 4 * params->n_pheno + params->nChrom + params->ncov;

a1 = 4.0 * num_phenotype + num_chromosome + num_covariate

# total_ram += std::max( params->block_size + params->n_pheno * params->n_ridge_l0 * t_eff, p_eff * params->n_ridge_l0 * b_eff );

b1 = block_size + num_phenotype * num_ridge_l0 * t_eff
b2 = p_eff * num_ridge_l0 * b_eff
max_b = max(b1, b2)

# total_ram *= params->n_samples * sizeof(double);

sample_cost = (a1 + max_b) * num_sample * 8

cat("The momory for storing sample level data:", ceiling(sample_cost / 1000^3), "GB\n")

# total_ram += params->nvs_stored * sizeof(struct snp);

variant_cost = num_variant * 300

cat("The momory for storing variant data:", ceiling(variant_cost / 1000^3), "GB\n")

# if( params->getCorMat ){ // M^2 (x2 with txt output) + 3NB
#     total_ram += (params->cor_out_txt && (params->ld_sparse_thr == 0) ? 2 : 1) * params->extract_vars_order.size() * params->extract_vars_order.size() * sizeof(double);
#     total_ram += params->n_samples * params->block_size * sizeof(double);
# }

other_cost = num_sample * block_size * 8

# if( params->use_loocv ) total_ram += params->chunk_mb * 1e6; // max amount of memory used for LOO computations involved

chunk_mb = 1000
other_cost += chunk_mb * 1000^2

# if( params->mask_loo ) total_ram += 1e9; // at most 1GB
other_cost += 1000^3

cat("The momory for storing other data:", ceiling(other_cost / 1000^3), "GB\n")

# if( params->vc_test ) total_ram += 2 * params->max_bsize * params->max_bsize * sizeof(double); // MxM matrices
# total_ram /= 1000.0 * 1000.0; 
# if( total_ram > 1000 ) {
#   total_ram /= 1000.0; 
#   ram_unit = "GB";
# } else ram_unit = "MB";

final_memory_gb = ceiling((sample_cost + variant_cost + other_cost) / 1000^3 * 1.3)

cat("The final memory with factor 1.3 cost:", final_memory_gb, "GB\n")

writeLines(paste0(final_memory_gb), "final_memory_gb.txt")

EOF

R -f script.r

  >>>

  runtime {
    cpu: 1
    docker: "~{docker}"
    preemptible: 1
    disks: "local-disk 5 HDD"
    memory: "1 GiB"
  }

  output {
    Int step1_memory_gb = read_int("final_memory_gb.txt")
  }
}

task Regenie4Step1FitModel {
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

    String step1_option = "--loocv --bsize 1000"

    String output_prefix

    Int memory_gb = 20
    Int? memory_gb_override
    Int cpu = 10

    Float disk_size_factor=2

    #String docker = "skoyamamd/regenie:3.4.2"
    #String docker = "quay.io/biocontainers/regenie:4.0--h90dfdf2_1"
    String docker = "shengqh/regenie4:20241127"
  }

  Int disk_size = ceil(size([input_pgen, input_pvar, input_psam], "GB") * disk_size_factor) + 10

  Int final_memory_gb = select_first([memory_gb_override, memory_gb])

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

mv ~{output_prefix}_pred.list old.list
 
while IFS=' ' read -r phenotype old_loco_file; do
  new_loco_file="~{output_prefix}.${phenotype}.loco"
  echo mv ${old_loco_file} ${new_loco_file}
  mv ${old_loco_file} ${new_loco_file}
  echo ${phenotype} ${new_loco_file} >> ~{output_prefix}_pred.list
done < old.list

>>>

  runtime {
    docker: docker
    preemptible: 1
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: final_memory_gb + " GiB"
  }
  output {
    File pred_list_file = "~{output_prefix}_pred.list" 
    Array[File] pred_loco_files = glob("*.loco")
  }
}

task Regenie4Step2AssociationTest {
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

    String? chromosome

    String output_prefix

    Int memory_gb = 20
    Int cpu = 10

    #String docker = "skoyamamd/regenie:3.4.2"
    #String docker = "quay.io/biocontainers/regenie:4.0--h90dfdf2_1"
    String docker = "shengqh/regenie4:20241127"
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

# Modified from https://github.com/briansha/Regenie_WDL/blob/master/regenie.wdl
task MergeRegenieChromosomeResults {
  input {
    Array[File] regenie_chromosome_files
    Array[String] phenotype_names
    Array[String] chromosome_list

    String regenie_prefix

    String output_prefix

    String docker = "ubuntu:20.04"
    Float memory = 3.5
    Int? disk_size_override
    Int cpu = 1
    Int preemptible = 1
    Int maxRetries = 0
  }
  Float regenie_files_size = size(regenie_chromosome_files, "GiB")
  Int disk = select_first([disk_size_override, ceil(30.0 + 2.0 * regenie_files_size)])

  String chr1 = chromosome_list[0]
  String pheno1 = phenotype_names[0]
  String regenie_file1 = "~{regenie_prefix}.chr~{chr1}_~{pheno1}.regenie"

command <<<

set -euo pipefail

for file in ~{sep=' ' regenie_chromosome_files}; do 
  mv $file .
done

for pheno in ~{sep=' ' phenotype_names}; do 
  head -n 1 ~{regenie_file1} > ~{output_prefix}.$pheno.regenie
  for chr in ~{sep= ' ' chromosome_list}; do 
    tail -n +2 ~{regenie_prefix}.chr${chr}_${pheno}.regenie >> ~{output_prefix}.$pheno.regenie
    rm ~{regenie_prefix}.chr${chr}_${pheno}.regenie
  done 

  gzip ~{output_prefix}.$pheno.regenie
done

>>>

  runtime {
    docker: docker
    memory: memory + " GiB"
    disks: "local-disk " + disk + " HDD"
    cpu: cpu
    preemptible: preemptible
    maxRetries: maxRetries
  }
  output {
    Array[File] phenotype_regenie_files = glob("*.regenie.gz")
  }
}

# Modified from https://github.com/briansha/Regenie_WDL/blob/master/regenie.wdl
# QQ and Manhattan plots
task RegeniePlots {
  input {
    File regenie_file
    String output_prefix

    # Runtime
    String docker = "shengqh/report:20241126"
    Float memory = 16.0
    Int? disk_size_override
    Int cpu = 1
    Int preemptible = 1
    Int maxRetries = 0
  }
  Float regenie_files_size = size(regenie_file, "GiB")
  Int disk_size = select_first([disk_size_override, ceil(10.0 + regenie_files_size)])

  # Plots are produced for each phenotype.
  # For each phenotype, a file containing all of the hits from Step 2 is output.
  # For each phenotype, a file containing a subset of all of the hits where "-LOG10P > 1.3" from Step 2 is output.
  command <<<
set -euo pipefail

cat <<EOF > script.r

library(data.table)
library(qqman)

regenie_output <- fread("~{regenie_file}")
regenie_ADD_subset <-subset.data.frame(regenie_output, TEST=="ADD")
regenie_ADD_subset[,"CHROM"] <-as.numeric(unlist(regenie_ADD_subset[,"CHROM"]))
regenie_ADD_subset[,"LOG10P"] <-as.numeric(unlist(regenie_ADD_subset[,"LOG10P"]))
regenie_ADD_subset[,"GENPOS"] <-as.numeric(unlist(regenie_ADD_subset[,"GENPOS"]))

png("~{output_prefix}.qqplot.png", width=5, height=5, units="in", res=300)
p = 10 ^ (-1 * (as.numeric(unlist(regenie_ADD_subset[,"LOG10P"]))))
print(qq(p))
dev.off()

png("~{output_prefix}.manhattan.png", width=10, height=5, units="in", res=300)
print(manhattan(regenie_ADD_subset, chr="CHROM", bp="GENPOS", snp="ID", p="LOG10P", logp=FALSE, annotatePval = 1E-5))
dev.off()

EOF

R -f script.r

>>>
  runtime {
    docker: docker
    memory: memory + " GiB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    preemptible: preemptible
    maxRetries: maxRetries
  }
  output {
    File qqplot_png = "~{output_prefix}.qqplot.png"
    File manhattan_png = "~{output_prefix}.manhattan.png"
  }
}