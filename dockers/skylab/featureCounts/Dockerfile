FROM python:3.6.2

LABEL maintainer="Lantern Team <lantern@broadinstitute.org>" \
      software="subread package" \
      version="2.0.1" \
      description="RNA-seq high-performance read alignment, quantification and mutation discovery" \
      website="http://subread.sourceforge.net/"

# Install compiler 
RUN apt-get update --fix-missing && apt-get install -y wget 

COPY requirements.txt .
RUN pip3 install -r requirements.txt    

# Install subread 
WORKDIR /usr/local/ 
ENV VERSION="2.0.1"
RUN wget "https://downloads.sourceforge.net/project/subread/subread-${VERSION}/subread-${VERSION}-source.tar.gz" \
     && tar -xzvf subread-${VERSION}-source.tar.gz
WORKDIR /usr/local/subread-${VERSION}-source/src 
RUN make -f Makefile.Linux 
ENV PATH /usr/local/subread-${VERSION}-source/bin/:$PATH
# Cleanup
RUN apt-get clean

# copy the script that removes alignments spanning intron-exon junctions
RUN mkdir /tools
WORKDIR /tools
COPY remove-reads-on-junctions.py .
