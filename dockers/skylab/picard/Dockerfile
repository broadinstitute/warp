FROM openjdk:8-jre

LABEL maintainer="Green Lantern <lantern@broadinstitute.org>" \
  software="Picard" \
  version="2.11.1" \
  description="A set of command line tools (in Java) for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF." \
  website="http://broadinstitute.github.io/picard"

# Please follow the below instructions to invoke picard when you are using this docker image:
# java jvm-args -jar /usr/picard/picard.jar PicardToolName OPTION1=value1 OPTION2=value2...
ENV version 2.20.4
WORKDIR /usr/picard
ADD https://github.com/broadinstitute/picard/releases/download/${version}/picard.jar ./picard.jar

# Install R dependencies for Picard
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y r-base

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('reshape')"
RUN Rscript -e "install.packages('gplots')"
RUN Rscript -e "install.packages('ggplot2')"
