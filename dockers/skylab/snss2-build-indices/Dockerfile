FROM python:3.6.2

LABEL maintainer="Farzaneh Khajouei <fkhajoue@broadinstitute.org>" \
      software="subread package" \
      version="2.0.2" \
      description="RNA-seq high-performance read alignment, quantification and mutation discovery" \
      website="http://subread.sourceforge.net/"

# Install compiler 
RUN apt-get update --fix-missing && apt-get install -y \
  wget 

COPY requirements.txt .
RUN pip3 install -r requirements.txt    

# Install subread 
WORKDIR /usr/local/ 

# copy the script that removes alignments spanning intron-exon junctions
RUN wget https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz && \
  tar -xf 2.7.8a.tar.gz
RUN chmod +x /usr/local/STAR-2.7.8a/bin/Linux_x86_64_static/STAR
ENV PATH /usr/local/STAR-2.7.8a/bin/Linux_x86_64_static/:$PATH
RUN mkdir /script
WORKDIR /script


COPY add-introns-to-gtf.py .
COPY modify_gtf_human.sh .
RUN chmod +x modify_gtf_human.sh
COPY modify_gtf_mouse.sh .
RUN chmod +x modify_gtf_mouse.sh
ENV PATH /script/:$PATH
