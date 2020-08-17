#!/bin/bash

docker build -t quay.io/humancellatlas/snap-breakout:0.0.1 .
docker push quay.io/humancellatlas/snap-breakout:0.0.1
