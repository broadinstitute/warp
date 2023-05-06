cd /nobackup/h_cqs/shengq2/test_data

gatk SplitSamByNumberOfReads --INPUT P100.0.unmapped.bam --OUTPUT . -N_READS 20000000 -TOTAL_READS 100245506 --OUT_PREFIX P100.0.unmapped

total_reads=$(samtools view -c P100.0.unmapped_0001.bam)

gatk SplitSamByNumberOfReads --INPUT P100.0.unmapped_0001.bam --OUTPUT . -N_READS 2000000 -TOTAL_READS 16707586 --OUT_PREFIX P100.0.unmapped_0001

gatk SplitSamByNumberOfReads --INPUT P100.0.unmapped_0001_0001.bam --OUTPUT . -N_READS 700000 -TOTAL_READS 1856400 --OUT_PREFIX P100.0.unmapped_0001_0001