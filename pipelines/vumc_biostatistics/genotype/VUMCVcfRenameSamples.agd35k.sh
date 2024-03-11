cd /nobackup/h_cqs/shengq2/biovu/cromwell

# run workflow with correct header
java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCVcfRenameSamples.wdl \
  -i /home/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCVcfRenameSamples.agd35k.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json

