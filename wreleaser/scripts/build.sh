#!/usr/bin/env bash

platforms=("linux/amd64" "darwin/amd64")

for platform in "${platforms[@]}"
do
    platform_split=($(echo $platform | tr "/" "\n"))
    GOOS=${platform_split[0]}
    GOARCH=${platform_split[1]}
    output_name='wreleaser-'$GOOS'-'$GOARCH

    env GOOS=$GOOS GOARCH=$GOARCH go build -o $output_name ../
    if [ $? -ne 0 ]; then
        echo 'An error has occurred! Aborting the script execution...'
        exit 1
    fi
done
