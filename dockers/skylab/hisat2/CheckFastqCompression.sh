#!/bin/bash
# fix names if necessary.
fastq="$1"
echo $fastq

if (file $fastq | grep -q compressed); then
    if [[ $fastq != *.gz ]]; then
        if [[ $fastq != *.fastq ]]; then
            FQ=$fastq.fastq.gz
            mv $fastq $fastq.fastq.gz
        else
            FQ=$fastq.gz
            mv $fastq $fastq.gz
        fi
    else
        FQ=$fastq
    fi
fi

echo $FQ
