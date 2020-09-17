# 1.5.0
2020-08-18

* Added a meta section with 'allowNestedInputs' set to 'true' to allow the workflow to use with nested inputs with Cromwell 52

# 1.4
2020-06-23

* Updated to use new version of CombineGenotypingArrayVcfs to remove AutoCallCallRate from 
single sample VCF before merging

# 1.3

2020-04-08

* Fixes arrays data delivery reference skew due to references moving to google buckets

# 1.2

2020-02-25

* Update wdl to use updated CombineGenotypingArrayVcfs (migrated to picard (public) repo).

# 1.1

2019-11-21

* Renamed file 'MultiSampleArraysWf.wdl' to 'MultiSampleArrays.wdl'
* Added a specified docker for each task in MultiSampleArrays

# 1.0
Initial release of the Multi-sample Arrays pipeline
