cd /nobackup/h_cqs/shengq2/biovu/cromwell

# run workflow with correct header
java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /nobackup/h_cqs/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCVcfExtractSamples.wdl \
  -i /nobackup/h_cqs/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCVcfExtractSamples.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json

