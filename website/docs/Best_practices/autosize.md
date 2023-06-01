---
sidebar_position: 2
---

# Autosizing disk for Google Cloud
When writing a WDL workflow, you have the option to customize your runtime environments, like the size of a virtual computer's disk, memory, and the type of disk you want to use. While you can hardcode standard values for these runtime parameters, there are times when dynamically calculating your parameters based on file sizes could be useful, like when you’re working with large data sets. If you find yourself having to modify your runtime values often, it’s probably beneficial to use some autosizing features.

Below, the Broad Pipeline Development team shares some tips and tricks for autosizing disk size for Google Cloud VMs.

## Autosizing disk based on file size
One simple step for autosizing the disk is to dynamically calculate a target disk size based on the size of your different input files. There are some useful WDL functions for this process, like `size()` and `ceil()`. The `size()` function returns the size of a designated file, whereas the `ceil()` function rounds up to the nearest integer. Combining these functions allows you to obtain an integer that reflects the file’s size; you can then assign that output to a variable that can be used throughout the different WDL tasks. 

An example of using autosize functions is shown in the code below; we combine the size() and ceil() functions to return the size of an input file called `input_file`:

```wdl
Int dynamic_disk_size = ceil(size(input_file,"GiB")) + 20
```
The `GiB` specified after the input file is just a size unit. You can see the full list of available size units in the [WDL 1.0 spec’s float size description](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#float-sizefile-string). The constant added to the end helps with small files sizes (for example, asking for only 1 GB of disk can be a problem for Google Cloud VMs).

The resulting output file size is assigned to the integer `dynamic_disk_size`. This variable can then be used to specify the disk size in the WDL’s runtime section, as shown below:   

```wdl
runtime {
  disks: "local-disk " + dynamic_disk_size + " HDD"
}
```

## Accounting for extra disk space and intermediate files 
At times, specifying disk size using an input file’s size alone may not meet the disk requirements. For example, if your workflow generates intermediate files, the disk needs to be large enough to handle these as well. You can account for intermediate files by setting a multiplier ( an integer) to an input variable. That way you can multiply the dynamically calculated file size by the number of intermediate files you generate, assuming they are roughly the same size as the input. The code below demonstrates this with the input `disk_multiplier`, which is set to 2 to handle a single intermediate file:

```wdl
 input {
  File clinical_input
  Int disk_multiplier = 2
}
  Int dynamic_disk_size = ceil(size(clinical_input,"GiB"))*disk_multuplier  + 20

runtime {
  disks: "local-disk " + dynamic_disk_size + " HDD"
}
```
A good rule of thumb for maximizing disk efficiency is to only use up to 80% of your disk. If you need to add extra space even after accounting for input and intermediate files, you can set up an input variable for extra disk space. The code below uses the input `extra_disk` to add 500 GiB to the disk size dynamically calculated by the input file’s size.

```wdl
input {
  File clinical_input
  Int extra_disk = 500
}
  Int dynamic_disk_size = ceil(size(clinical_input,"GiB"))*2 + extra_disk
```

## Making optional disk input
You may want to add some flexibility to your workflow so that you can manually specify disk size in addition to dynamically calculating it. In this case, you’ll want to use an optional input and take advantage of the WDL function `select_first()`. This function selects the first defined function and returns it. 

Let’s take a look at the example code below, which sets up an optional input, `disk_size_gb`.  The code uses the WDL 1.0 `size()` and `ceil()` functions to assign the size of the file `clinical_input ` to a dynamically calculated input variable called `dynamic_disk_size`. It then uses `select_first()` to choose either the optional `disk_size_gb` value (if defined) or the dynamically calculated value.

```wdl
input {
  File clinical_input
  Int? disk_size_gb
}
  Int dynamic_disk_size = ceil(size(clinical_input,"GiB"))*2 + 500
  Int disk_size = select_first([disk_size_gb, dynamic_disk_size])
runtime {
  disks: "local-disk " + disk_size + " HDD"
}
```
The code above saves the output of the `select_first()` function to the variable `disk_size`, which is then used to assign the disk size in the WDL’s runtime section.

## Additional Considerations

### Disk types
The disk types you select can vary in cost and therefore affect how much wiggle room you want to add for your disk and memory sizes. The different disk types are described in [Cromwell documentation](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#disks), but basically, you can select between a more expensive solid-state drive (SSD) or a standard hard disk drive HDD. SSD is useful for efficient input/output streaming. If you opt for an SSD disk, you’ll want to restrict the amount of disk size padding and multipliers you use. Because SSD can work more efficiently for some tasks, you might not require the 500 GB extra that you include when working with HDD disks. 

### Google Cloud vs. other platforms
The suggestions in this article are specific to Google Cloud VMs; other platforms may only have disk_size as an option or different types of disks. Remember that these suggestions may not apply across different platforms.

## Resources
* [Ceil function description](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#int-floorfloat-int-ceilfloat-and-int-roundfloat)
* [Size function description](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#float-sizefile-string)
* [Select_first function description](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#x-select_firstarrayx)
* [Runtime defaults defined by Cromwell](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#runtime-attribute-descriptions)
