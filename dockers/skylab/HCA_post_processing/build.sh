#!/bin/bash

tag=$1

if [ -z $tag ]; then
    echo -e "\nYou must provide a tag"
    echo -e "\nUsage: bash build_docker.sh TAG\n"
    exit 1
fi

docker build . --tag=quay.io/humancellatlas/hca_post_processing:$tag

echo "You can now push with docker push quay.io/humancellatlas/hca_post_processing:$tag"
