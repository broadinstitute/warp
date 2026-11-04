---
sidebar_position: 6
---

# Task execution - tips for using the WDL task command section
Every WDL task has a [command section](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#command-section) where you can call software tools and specify parameters to help transform your data into meaningful output. This section is like a terminal for whatever environment you’re using to execute your WDL script. That environment can be a virtual computer set up by a Docker container or it can be your local computer. If you’re using Cromwell to execute your WDL (as happens in the cloud-based platform [Terra](https://app.terra.bio/)), the command section is run after Cromwell has resolved the task inputs but before it assesses outputs. 

The actual commands you use in the task command section depend on what tools and operating systems are available in the execution environment. [WARP](https://github.com/broadinstitute/warp/tree/master) pipelines, for example, often set up virtual machines using Docker containers with Alpine Linux-based operating systems. This means that the command section should contain commands that work in Alpine. If additional software is installed on top of Alpine, that software's commands will also work. WARP workflows often require custom python scripts and that’s why python is installed on top of WARP’s Alpine-based (or other OS) Dockers. Python is one example, but you can install any language on a docker and then use language-specific commands or scripts from the WDL task command section. In addition to these commands, you can also point to paths for software (such as the path to a jar in the Docker) as well as input/output files that are in the Docker container.

If a script is small, you can even write scripts in line. The [Resources](#resources) section below links to examples of workflow code that use in-line scripting. If you’re using in-line scripting, remember to keep your code tidy (check out the WARP WDL formatting suggestions). And remember, commands will only work if installed on the environment in which the WDL is running; what works on your local machine will not necessarily work in a virtual machine built off a Docker and vice-versa.  

If you read the WDL 1.0 spec’s [Command Section](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#command-section), you’ll notice there are two options for specifying a command section: curly brackets `{...}` or angled brackets `<<<...>>>`. This choice impacts how you specify inputs from your WDL inputs section. The curly brackets use `${...}` or `~{...}` to specify a WDL input whereas the angled brackets explicitly use `~{...}`. If you’re familiar with bash, you might recognize a problem with using the dollar sign to specify WDL inputs; bash also uses a dollar sign for inputs. If you specify inputs using `${...}`, it will interpret the input as a WDL input and not a bash variable. For this reason, the Broad Pipeline Development team recommends using the angled brackets option, allowing you to avoid collisions and specify bash variables with dollar signs.

## Resources
Below are some example WDLs that either use scripting languages in-line or call scripts written in particular languages.

* [Example of WDL that uses in-line perl script (see line #348 of the Build Indices workflow)](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/build_indices/BuildIndices.wdl)
* [Example of WDL task that uses in-line python script (see line #28 of Utilities WDL](https://github.com/broadinstitute/warp/blob/master/tasks/wdl/Utilities.wdl)
* [Example of WDL that uses an R script (line 39 of the emptyDrops task)](https://github.com/broadinstitute/warp/blob/master/tasks/wdl/RunEmptyDrops.wdl)
