FROM gcr.io/google.com/cloudsdktool/cloud-sdk

LABEL maintainer="Lantern Team <lantern@broadinstitute.org>"

RUN pip3 install --upgrade pip

COPY requirements.txt .
RUN pip3 install numpy==1.17.0
RUN pip3 install cython==0.29.15
RUN pip3 install -r requirements.txt

RUN mkdir /tools
WORKDIR /tools

COPY HCA_create_adapter_json.py .
COPY file_utils.sh .
COPY create_input_metadata_json.py .
