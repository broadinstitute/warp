# Build docker image from Dockerfile in this directory
# This is to help keep track of versions, etc.
# Wraps around ImputationPipeline/build_push_docker.sh.

dockerfile_directory=imputation_bcftools_vcftools_docker
image_version=v1.0.0 # as of Jan 25 2021
ubuntu_version=20.04 # as of Jan 25 2021

wd=$(pwd)
cd "$(dirname $0)" || exit

../build_push_docker.sh \
  --directory ${dockerfile_directory} \
  --ubuntu-version ${ubuntu_version} \
  --image-version-tag ${image_version} \
#  --dry-run


cd "${wd}" || exit
