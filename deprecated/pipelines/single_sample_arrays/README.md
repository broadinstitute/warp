# Building the dockers
# Arrays uses three docker images

# arrays-picard
The arrays-picard-private docker image contains the picard-private jar, vault, and gsutil.
It is built with the './build_arrays_picard_private_docker.sh' script.
The script builds it in a local (temp) docker directory
and then pushed to us.gcr.io/broad-arrays-prod/arrays-picard-private

# illumina-autocall
The illumina-autocall docker image contains Illumina's autocall tool
It is built with the './build_illumina_autocall_docker.sh' script.
The script builds it in a local (temp) docker directory
and then pushed to us.gcr.io/broad-arrays-prod/illumina-autocall

# zcall
The zcall docker image contains the python script Zcall
It is built with the './build_zcall_docker.sh' script.
The script builds it in a local (temp) docker directory
and then pushed to us.gcr.io/broad-arrays-prod/zcall

# Other
When building these docker images you may be asked for your unix password
for the `scp`s used when building the images.

