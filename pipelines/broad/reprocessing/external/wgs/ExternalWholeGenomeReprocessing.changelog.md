# 1.1.1
2020-10-01

* Removed extra trailing slash in ouput directory from cloud to cloud copy job

# 1.1.0
2020-08-18

* Added a meta section with 'allowNestedInputs' set to 'true' to allow the workflow to use with nested inputs with Cromwell 52

# 1.0.2
2020-07-31

* Update various tasks to stop using phusion/baseimage:latest docker image (it has been removed).  Start using a Google-hosted base image in it's stead.

# 1.0.1
2020-07-15

* Remove GetBWAVersion as a task and moved it to SamToFastqAndBwaMemAndMba

# 1.0
2020-04-30

* Initial release of the ExternalWholeGenomeReprocessing pipeline
