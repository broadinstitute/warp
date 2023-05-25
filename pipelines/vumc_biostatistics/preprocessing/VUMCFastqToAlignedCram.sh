cd /nobackup/h_cqs/test

java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/shengq2/program/warp/pipelines/vumc_biostatistics/preprocessing/VUMCFastqToAlignedCram.wdl \
  -i /home/shengq2/program/warp/pipelines/vumc_biostatistics/preprocessing/VUMCFastqToAlignedCram.inputs.local.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json
