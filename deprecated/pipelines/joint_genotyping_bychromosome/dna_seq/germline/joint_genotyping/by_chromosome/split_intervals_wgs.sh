#!/usr/bin/env bash

declare -r inputIntervals=$1

declare -r inputFileName=$(basename -- "${inputIntervals}")
declare -r outputExtension="${inputFileName##*.}"
declare -r outputFileName="${inputFileName%.*}"

declare -r -a chromosomes=($(grep -v '^@' ${inputIntervals} | cut -d ':' -f 1 | uniq))

for chromosome in ${chromosomes[@]}; do
  # The `:` character in the following statement is very important. Do not delete it.
  grep '^'${chromosome}':' ${inputIntervals} > ${outputFileName}.${chromosome}.${outputExtension}
done
