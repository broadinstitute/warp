cd /nobackup/h_vangard_1/shengq2/biovu

if [[ ! -s merged_chr2.vcf.gz ]]; then
  gsutil -m cp gs://original-plink/working_set_redeposit_megaex1-1/redeposit_vcf_imputed/TOPMed_imputed/remerged_EUR/merged_chr2.vcf.gz .
fi

java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCVcfIncludeSamples.wdl \
  -i /home/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCVcfIncludeSamples.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json