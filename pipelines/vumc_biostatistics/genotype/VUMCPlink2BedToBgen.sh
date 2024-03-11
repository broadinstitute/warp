if [ ! -d "/workspace/shengq2/temp/bgen" ]; then
  mkdir /workspace/shengq2/temp/bgen
fi

cd /workspace/shengq2/temp/bgen

java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCPlink2BedToBgen.wdl \
  -i /data/cqs/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCPlink2BedToBgen.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json