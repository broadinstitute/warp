cd /nobackup/h_cqs/shengq2/biovu/wgs

java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local.conf \
  -jar /data/cqs/softwares/wdl/cromwell-84.jar \
  run /home/shengq2/program/warp/pipelines/vumc_biostatistics/dna_seq/germline/single_sample/wgs/VUMCIlluminaCoverage.wdl \
  -i /home/shengq2/program/warp/pipelines/vumc_biostatistics/dna_seq/germline/single_sample/wgs/VUMCIlluminaCoverage.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json