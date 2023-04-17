cd /nobackup/h_cqs/test

java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /data/cqs/softwares/vumc_biostatistics/warp/pipelines/vumc_biostatistics/preprocessing/VUMCSamHeader.wdl \
  -i /data/cqs/softwares/vumc_biostatistics/warp/pipelines/vumc_biostatistics/preprocessing/VUMCSamHeader.inputs.local.unmapped.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json
