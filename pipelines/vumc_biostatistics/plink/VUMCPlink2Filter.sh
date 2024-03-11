cd /nobackup/h_cqs/shengq2/test
java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /nobackup/h_cqs/shengq2/program/warp/pipelines/vumc_biostatistics/plink/VUMCPlink2Filter.wdl \
  -i /nobackup/h_cqs/shengq2/program/warp/pipelines/vumc_biostatistics/plink/VUMCPlink2Filter.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json