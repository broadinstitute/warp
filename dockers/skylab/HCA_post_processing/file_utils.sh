#!/bin/bash

function get_timestamp() {
  local -r input_file=${1}
  timestamp=$(gsutil ls -l ${input_file} | egrep -o "([0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z)")
  echo ${timestamp}
}

function get_crc() {
  local -r input_file=${1}
  crc=$(gsutil hash -h ${input_file} | awk '/crc32c/ { print $3 }')
  echo ${crc}
}

function get_size() {
  local -r input_file=${1}
  size=$(gsutil stat ${input_file} | awk '/Content-Length/ { print $2 }')
  echo ${size}
}
