cd /nobackup/h_cqs/shengq2/biovu/agd35k/ICA-AGD/small_callsets

java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /nobackup/h_cqs/shengq2/program/warp/pipelines/vumc_biostatistics/plink/VUMCPlink2FilterRegion.wdl \
  -i /nobackup/h_cqs/shengq2/program/warp/pipelines/vumc_biostatistics/plink/VUMCPlink2FilterRegion.exome.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json

java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /nobackup/h_cqs/shengq2/program/warp/pipelines/vumc_biostatistics/plink/VUMCPlink2FilterRegion.wdl \
  -i /nobackup/h_cqs/shengq2/program/warp/pipelines/vumc_biostatistics/plink/VUMCPlink2FilterRegion.clinvar.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json  