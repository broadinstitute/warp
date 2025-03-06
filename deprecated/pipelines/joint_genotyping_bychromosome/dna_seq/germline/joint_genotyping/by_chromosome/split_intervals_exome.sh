#!/usr/bin/env bash

declare -r inputIntervalList=$1

declare -r inputFileName=$(basename -- "${inputIntervalList}")
declare -r outputExtension="${inputFileName##*.}"
declare -r outputFileName="${inputFileName%.*}"

declare -r -a chromosomes=($(grep -v '^@' ${inputIntervalList} | cut -f 1 | uniq))

declare -r headerTmp=$(mktemp)

grep '^@' ${inputIntervalList} > ${headerTmp}

for chromosome in ${chromosomes[@]}; do
  intervalTmp=$(mktemp)
  grep '^'${chromosome}'\t' ${inputIntervalList} > ${intervalTmp}
  cat ${headerTmp} ${intervalTmp} > ${outputFileName}.${chromosome}.${outputExtension}
  rm ${intervalTmp}
done

rm ${headerTmp}

