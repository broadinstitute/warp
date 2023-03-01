cd /scratch/cqs/shengq2/test

java -Dconfig.file=/home/shengq2/program/cqsperl/config/wdl/cromwell.local.conf -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/shengq2/program/warp/pipelines/vumc_biostatistics/preprocessing/VUMCFastqQC.wdl \
  -i /home/shengq2/program/warp/pipelines/vumc_biostatistics/preprocessing/VUMCFastqQC.inputs.local.json \
  --options ~/program_other/si2022/cromwell.options.json
