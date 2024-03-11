cd /nobackup/h_vangard_1/shengq2/biovu

if [[ ! -s EUR_chr22.n100.vcf.gz ]]; then
  zcat /nobackup/h_cqs/shengq2/biovu/megaex1.1/TOPMed_imputed/EUR/merged_chr22.vcf.gz | head -n 100 > EUR_chr22.n100.vcf
  bgzip EUR_chr22.n100.vcf
  tabix -p vcf EUR_chr22.n100.vcf.gz
fi

# run workflow with wrong header
java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCVcfIncludeSamples.wdl \
  -i /home/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCVcfIncludeSamples.header_wrong.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json
