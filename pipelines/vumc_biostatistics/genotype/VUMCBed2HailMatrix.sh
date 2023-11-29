if [ ! -d "/workspace/shengq2/temp/hail" ]; then
  mkdir /workspace/shengq2/temp/hail
fi

cd /workspace/shengq2/temp/hail

java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCBed2HailMatrix.wdl \
  -i /data/cqs/shengq2/program/warp/pipelines/vumc_biostatistics/genotype/VUMCBed2HailMatrix.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json