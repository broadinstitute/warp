cd /scratch/cqs/jstolze/terra/CountCram_testing/

java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/stolzelk/warp/pipelines/vumc_biostatistics/preprocessing/VUMCCramQC2.wdl \
  -i /home/stolzelk/warp/pipelines/vumc_biostatistics/preprocessing/VUMCCramQC2.inputs.local.mapped.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json