cd /scratch/cqs/shengq2/temp/

java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/shengq2/program/warp/pipelines/vumc_biostatistics/preprocessing/VUMCCramQC2.wdl \
  -i /home/shengq2/program/warp/pipelines/vumc_biostatistics/preprocessing/VUMCCramQC2.inputs.local.mapped.small.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json