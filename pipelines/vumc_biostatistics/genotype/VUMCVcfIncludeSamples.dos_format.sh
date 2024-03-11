cd /nobackup/h_vangard_1/shengq2/biovu

if [[ ! -s EAS_chr20.n100.vcf.gz ]]; then
  zcat /nobackup/h_cqs/shengq2/biovu/megaex1.1/TOPMed_imputed/EAS/chr20.dose.vcf.gz | head -n 100 > EAS_chr20.n100.vcf
  bgzip EAS_chr20.n100.vcf
  tabix -p vcf EAS_chr20.n100.vcf.gz
fi

# run workflow with correct header
java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCVcfIncludeSamples.wdl \
  -i /home/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCVcfIncludeSamples.dos_format.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json

