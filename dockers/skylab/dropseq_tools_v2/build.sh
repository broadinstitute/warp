#!/bin/bash

tag=2.3.0
image="quay.io/humancellatlas/secondary-analysis-dropseqtools"

if [ -z $tag ]; then
    echo -e "\nYou must provide a tag"
    echo -e "\nUsage: bash build_docker.sh TAG\n"
    exit 1
fi

docker build -t $image:$tag .

echo "You can now push with docker push $image:$tag"
