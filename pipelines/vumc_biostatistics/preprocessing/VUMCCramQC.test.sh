cd /scratch/cqs/test

java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /data/cqs/softwares/vumc_biostatistics/warp/pipelines/vumc_biostatistics/preprocessing/VUMCCramQC.wdl \
  -i /data/cqs/softwares/vumc_biostatistics/warp/pipelines/vumc_biostatistics/preprocessing/VUMCCramQC.inputs.local.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json
