## Use hail as the base image, this is a clone
## We do not want to hit dockerhub because of the pull limit
## So we don't do anything here
FROM hailgenetics/hail:0.2.134-py3.11
COPY scripts/* . 
CMD ["bash"]
