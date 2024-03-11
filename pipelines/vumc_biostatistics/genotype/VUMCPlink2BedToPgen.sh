if [ ! -d "/nobackup/h_cqs/shengq2/biovu/agd35k/pgen" ]; then
  mkdir /nobackup/h_cqs/shengq2/biovu/agd35k/pgen
fi

cd /nobackup/h_cqs/shengq2/biovu/agd35k/pgen

java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCPlink2BedToPgen.wdl \
  -i /data/cqs/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCPlink2BedToPgen.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json