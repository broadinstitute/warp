mkdir -p /nobackup/h_cqs/shengq2/biovu/agd35k/pgen
cd /nobackup/h_cqs/shengq2/biovu/agd35k/pgen
java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /nobackup/h_cqs/shengq2/program/warp/pipelines/vumc_biostatistics/GWAS/VUMCQCFilterAndMergePgen.wdl \
  -i /nobackup/h_cqs/shengq2/program/warp/pipelines/vumc_biostatistics/GWAS/VUMCQCFilterAndMergePgen.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json