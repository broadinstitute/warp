cd /nobackup/h_cqs/shengq2/biovu/megaex/EUR_NewImputation

if [[ ! -s EUR-chr22.n100.vcf.gz ]]; then
  zcat EUR-chr22.dose.mac11.R203.vcf.gz | head -n 100 > EUR-chr22.n100.vcf
  bgzip EUR-chr22.n100.vcf
  tabix -p vcf EUR-chr22.n100.vcf.gz
fi

# run workflow with correct header
java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCVcfIncludeSamples.wdl \
  -i /nobackup/h_cqs/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCVcfIncludeSamples.NewEUR.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json

