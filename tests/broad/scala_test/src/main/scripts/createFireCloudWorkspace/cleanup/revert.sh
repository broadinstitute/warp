#!/usr/bin/env bash

gsutil -m cp -r gs://broad-gotc-prod-storage/pipeline/G96830/NA12878/v470/* gs://broad-gotc-dev-storage/pipeline/G96830/NA12878/v470/
gsutil -m cp -r gs://broad-gotc-prod-storage/pipeline/G94982/NA12891/v9/* gs://broad-gotc-dev-storage/pipeline/G94982/NA12891/v9/
gsutil -m cp -r gs://broad-gotc-prod-storage/pipeline/G94982/NA12892/v8/* gs://broad-gotc-dev-storage/pipeline/G94982/NA12892/v8/

java -Dclio.server.hostname=clio.gotc-dev.broadinstitute.org -Dclio.server.port=443 -Dclio.server.use-https=true -jar /seq/software/clio-client/current/clio-client.jar add-cram --data-type WGS --location GCP --project G96830 --sample-alias NA12878 --version 470 --metadata-location G96830_NA12878_470.json

java -Dclio.server.hostname=clio.gotc-dev.broadinstitute.org -Dclio.server.port=443 -Dclio.server.use-https=true -jar /seq/software/clio-client/current/clio-client.jar add-cram --data-type WGS --location GCP --project G94982 --sample-alias NA12891 --version 9 --metadata-location G94982_NA12891_9.json

java -Dclio.server.hostname=clio.gotc-dev.broadinstitute.org -Dclio.server.port=443 -Dclio.server.use-https=true -jar /seq/software/clio-client/current/clio-client.jar add-cram --data-type WGS --location GCP --project G94982 --sample-alias NA12892 --version 8 --metadata-location G94982_NA12892_8.json
