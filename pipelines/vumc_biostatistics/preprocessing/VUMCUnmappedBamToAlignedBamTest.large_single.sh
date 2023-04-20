cd /nobackup/h_cqs/shengq2/test/wdl

java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/shengq2/program/warp/pipelines/vumc_biostatistics/preprocessing/VUMCUnmappedBamToAlignedBam.wdl \
  -i /home/shengq2/program/warp/pipelines/vumc_biostatistics/preprocessing/VUMCUnmappedBamToAlignedBam.inputs.local.large_single.json \
  --options ~/program_other/si2022/cromwell.options.json
