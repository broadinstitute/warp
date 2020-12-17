#!/bin/bash
# A simple script to iterate over all wdl files in pipelines, projects, and beta-pipelines directories and run 'miniwdl check' on them.
# Do not exit on the first failure.

mainExitCode=0
for file in $(find pipelines/ projects/ -name '*.wdl' -type f -print)
do
  miniwdl check ${file}
  exitCode=$?
  if [[ ${exitCode} != 0 ]]; then
    mainExitCode=${exitCode}
  fi
done
exit ${mainExitCode}