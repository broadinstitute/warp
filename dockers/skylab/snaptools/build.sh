#!/bin/bash

docker build . -t quay.io/humancellatlas/snaptools:0.0.1

echo You can now push to quay.io using the following command
echo    docker push quay.io/humancellatlas/snaptools:0.0.1
