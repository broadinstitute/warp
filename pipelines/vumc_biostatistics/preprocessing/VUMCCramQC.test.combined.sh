cd /nobackup/h_cqs/jstolze/terra/CountCram_testing/incorperating_with_orig_QC/
java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/stolzelk/warp/pipelines/vumc_biostatistics/preprocessing/VUMCCramQC.wdl \
  -i /home/stolzelk/warp/pipelines/vumc_biostatistics/preprocessing/VUMCCramQC.inputs.local.mapped_combined.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json