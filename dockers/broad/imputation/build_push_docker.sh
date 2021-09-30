set -e

while test $# -gt 0; do
  case "$1" in
    -h|--help)
      echo "Build and push a docker image from specified Dockerfile."
      echo "Call from anywhere in the palantir-workflows repo. Supply the following options:"
      echo "Usage:                     ./build_push_docker.sh [OPTIONS]"
      echo "Options:"
      echo "    --directory/-d         (required) directory containing the Dockerfile (eg. imputation_eagle_docker)"
      echo "    --image-tag/-i         (required) tag for this image (suggest using the image version)"
      echo "    --ubuntu-version/-u    (optional) version of ubuntu base image (defaults to 20.04)"
      echo "    --no-push/-p           (optional) build but do not push"
      echo "    --dry-run/-r           (optional) dry run (no build or push; can use to inspect variables)"
      echo "    --no-cache/-c          (optional) build docker image with no cache"
      exit 1
      ;;
    -d|--directory)
      shift
      DOCKER_DIR=$1
      shift
      ;;
    -i|--image-version-tag)
      shift
      IMG_TAG=$1
      shift
      ;;
    -u|--ubuntu-version)
      shift
      UBUNTU=$1
      shift
      ;;
    -r|--dry-run)
      shift
      DRY="true"
      ;;
    -p|--no-push)
      shift
      PUSH="false"
      ;;
    -c|--no-cache)
      shift
      NOCACHE="--no-cache"
      ;;
    *)
      echo "Invalid argument. Use --help or -h flags for usage information."
      exit 1
      ;;
  esac
done


if [[ -z $DOCKER_DIR ]]; then
  echo "No docker path specified. Please specify a directory containing a dockerfile with -d or --directory."
  exit 1
fi

# check for missing arguments, fill defaults
if [[ -z $UBUNTU ]]; then
  echo "No ubuntu version specified. Using ${UBUNTU} as default."
  UBUNTU=20.04
fi
if [[ -z $IMG_TAG ]]; then
  echo "No image version specified. Please specify an image version with -i or --image-version."
  exit 1
fi
[[ -z "${DRY}" ]] && DRY=false
[[ -z "${PUSH}" ]] && PUSH=true


docker_dir=$(basename ${DOCKER_DIR})
docker_path=$(find .. -type d -name "${docker_dir}")
image_name=us.gcr.io/broad-dsde-methods/${docker_dir}:${IMG_TAG}
while true; do
  if [[ "${DRY}" == "true" ]]; then
    break;
  fi
  echo "This script will build and push ${image_name}. Do you want to proceed? (y/[n])"
  read yn
  [[ -z ${yn} ]] && yn=n
  case $yn in
    [Yy]* ) break;;
    [Nn]* ) exit 1;;
    * ) echo "Please answer yes or no.";;
  esac
done

# Execute commands only if this is not a dry run
function execute(){
        # Irrespective of whether dry run is enabled or not, we display
        # the command on the screen
    # shellcheck disable=SC2145
    echo "COMMAND: ${@}"
    # if dry run is enabled then simply return
    if [[ ${DRY} == "false" ]]; then
        eval "$@"
    fi
}

# Check the docker path


echo "Directory: ${docker_path}"
if ! [[ $(find "${docker_path}" -name Dockerfile) ]]; then
  echo "No Dockerfile found in this directory."
  exit 1
fi


image_name=us.gcr.io/broad-dsde-methods/${docker_dir}:${IMG_TAG}

echo "Ubuntu version: ${UBUNTU}"
echo "Image version tag: ${IMG_TAG}"



build_opts="-t ${image_name} --build-arg UBUNTU_VERSION=${UBUNTU} ${NOCACHE}"
execute "docker build ${docker_path} ${build_opts}"
execute "docker push ${image_name}"