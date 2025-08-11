---
sidebar_position: 1
---

# WDL cost optimization: Tips and tricks when working with Google Cloud in Terra

Reducing the cost of your WDL workflow is always a priority. Below, the Broad Institute’s Pipelines team provides some tips and tricks for optimizing workflow costs for runs using Google Cloud virtual machines (VM) from the bioinformatics platform, [Terra](https://app.terra.bio/). 

Overall, the majority of optimization comes down to understanding the size of your VM and how long you use it. Keep this in mind as you read the tips below and remember: one size does not necessarily fit all when it comes to optimizing your workflow.

## Tip 1: Know the difference between tools and bash commands

A big part of cost reduction is understanding what tools do the heavy lifting in your analysis because these tools will often dictate the size requirements of a VM and how long it runs. When we talk about software tools, we mean the time and memory-intensive tools required to analyze your data. For example, if you’re doing sequencing alignment, your tools might be bwa, STAR, or whatever your favorite aligner is. 

In contrast, bash commands are primarily used to help (re)format, (re)name, and move files around to get them ready to do the large analysis. These functions will likely have a smaller impact on the necessary VM size than your tools.

The distinction between tools and bash commands becomes important when we discuss running one tool per WDL task in Tip 2. Before reading on, think about what tools vs. bash operations you currently use for your analysis. 

## Tip 2: Modularize tools when possible

When you run a WDL workflow in the cloud, you specify a docker for each task, which will be used to set up a virtual machine (VM). Initiating a VM with a docker comes with some overhead costs so it might be tempting to avoid these by running all of your software tools in a single task, using only a single VM and docker. 

Unfortunately, the reality is that software tools aren’t 100% reliable, and when they fail, it means recreating a VM and rerunning your tools to successfully complete your workflow run. When you run multiple tools on a single VM or in a single task, it doesn’t matter if one tool fails or multiple tools fail, you’ll have to rerun everything that’s part of the WDL task. You’ll also have a hard time making updates to the code should you want to reuse one of the tools in a separate workflow later on.


For these reasons, we recommend modularizing workflows when possible, limiting the number of software tools you run per task. For example, if you're software tool takes more than 10 minutes to run, it's probably a good idea to put the tool in it's own task. When you modularize tools, you’ll be able to save your upstream analysis even if your tool fails, saving you the cost of time inside your VM. You’ll also be able to better customize the VM to meet the needs of the tools and run more efficiently.

In addition to modularizing tools, you can also consider whether you want to modularize bash commands. Bash commands are often less time intensive because they are quick steps that you use to prepare files for processing with your software tool. In the case of quick bash commands, it makes sense to group the, in the same WDL task as you would your tool as opposed to separating them into their own task. 

Let’s look at an example WDL task below where we first use bash commands to prepare sequencing data for genomic alignment with the STAR aligner, and then run the STAR tool:

```wdl
command <<<
  Set -e

# prepare reference with bash
mkdir genome_reference
Tar -xf “~{tar_star_reference}” -C genome_reference --strip component 1
Rm “~{tar_star_reference}

STAR \
-- runMode alignReads \
--runThreadN ~{cpu}
--genomeDir genome_reference
>>>
```

In the use case above, it makes sense to combine the bash operations and STAR. The bash operations are only used to make a new directory for the genome references and to compress files. After the files are ready, the heavy lifting of genome alignment is done by STAR.

Keep in mind there might be scenarios where splitting out more intensive bash operations might make sense. For example, if your bash operations are time intensive and your tool is not very reliable, you might want to save the outputs of your bash functions before moving on to running the tool. In this case, it would make sense to keep the bash commands and your tool in separate WDL tasks to avoid rerunning them both.

Similarly, there are exceptions to modularization of tools. For example, if you're working with a lot of intermediate files, it might make sense to string two tools together in a single task if it helps avoid costly egress. Always keep in mind timing and the reliability of your tools; remember that these tips are guidelines, not hard rules. 

## Tip 3: Avoid input/output timing costs (moving and loading lots of files)

Just like starting a VM has some overhead cost, localizing files also has some downsides, including using more time in your VM. And the more time you spend, the more your cost will go up. Your VM needs to find all your cloud files in order to run them through your different tools. That means if you have thousands of files, you’ll have to run initialization steps for each of them, which can take a bit of time and networking. 

Each time you move cloud files, you pay egress for the network transitions, so it’s important to find the balance in the number of files you decide to move. You have to weigh whether it costs more to move a large file vs. moving several smaller files. For example, you might find that it’s more cost-efficient to move a zipped \~100 GB file than to move 100, 1 GB files. After running a  test workflow, check your workflow logs to see what the timing is for localizing and moving files vs. running your tool. This might require trial and error when developing your workflow. 

## Tip 4: Run files in parallel when possible
If you need to run multiple files through a tool, you’ll have to decide whether to scatter those files across multiple VMs (running the tool in parallel), or run the files sequentially through the tool in one VM. 

Similar to the problem of running multiple tools per VM, running multiple files per VM also comes with the risk of rerunning your tool if a file should fail. If you’re running files that are large or prone to transient failures, it’s best to scatter them across VMs in parallel. 

Let’s take a look at an example WDL script that shows scattering; it performs variant calling on large BAM files with a tool called HaplotypeCaller. The WDL scatters the input array of BAM files so that the task runs on each individual BAM. It creates a single VM for each instance of the task. 

```wdl
version 1.0
workflow ScatterGatherExample {
  input {
    Array[File] sampleBAMs
  }

  scatter (sample in sampleBAMs) {
    call HaplotypeCallerERC { 
      input: 
        bamFile=sample 
    }
  }
  call GenotypeGVCF {
   input: 
     GVCFs = HaplotypeCallerERC.GVCF 
  }
}
```

If one of the BAM files in the array fails to run through HaplotypeCaller, the remaining files will still be stored. In consequence, you can rerun the WDL to analyze the failed file and not waste additional resources rerunning the files that already successfully ran through the task.

Another consideration when running files in parallel is to make sure that files are split into roughly equal size chunks. This is important if, for example, you have hard-coded values for memory and disk. In that case, each shard will get as much memory and disk as needed for your largest shard and the largest one can sometimes require a long time to run.

As with the other tips, running in parallel might not be necessary for all use cases. If you’re running smaller and fewer files, you can save money by running sequentially through the task in a single VM. 

## Tip 5: Use a VM that already has your reference
When you run WDL workflows in Terra, you have the option to use a VM that already has reference files loaded. This is important because you can avoid localization costs of necessary reference files, saving you time in the VM and ultimately, money.

Let’s look at an example of how using preloaded references helped optimize timing for the alignment task of a real workflow, the [Smart-seq2 Single Nucleus Multi-Sample Pipeline](https://github.com/broadinstitute/warp/blob/master/pipelines/wdl/smartseq2_single_nucleus_multisample/MultiSampleSmartSeq2SingleNucleus.wdl). 

Prior to using preloaded references, the timing for the alignment task using STAR was the following:

| Timing of task steps | Task step (from log file) |
| --- | --- |
| 00:52:29 | Starting container setup. |
| 00:52:38  | Localization script execution started. |
| 00:52:45  | Localizing input |
| 01:03:29 | Done localization. |
| 01:03:30 |  Running user action: docker |
| 01:17:14 |  ..... started STAR run |
| 01:17:15 |  ..... loading genome |
| 01:30:31 | ..... started mapping |
| 01:32:25 | ..... finished mapping |
| 01:32:27 |  ..... started sorting BAM |
| 01:32:29 | ..... finished successfully |
| 01:32:41 |  Starting delocalization. |
| 01:32:52 | Done delocalization. |

Notice that the localizing input takes ~ 11 min and that loading the genome takes an additional 13 min. In contrast, mapping the files to the genome (the actual task) takes only ~ 3 min. This means we spend most of our time just finding and loading files. 

In contrast, when we use preloaded references, the timing looks more like the following:

| Timing of task steps | Task step (from log file) |
| --- | --- |
| 00:52:29 | Starting container setup. |
| 00:52:38  | Localization script execution started. |
| 00:52:45  | Localizing input |
| 00:52:55 | Done localization. |

Notice it now only takes 10 seconds to localize the files! 

When you’re troubleshooting workflows in Terra, be sure to check out your workflow timing logs so you can identify which steps are costing you time in the VM.


## Conclusion
Workflow cost optimization is a bit of an art, but the payoff is worth the effort. Overall, remember that the size of your VM and the time you spend using it drive your cost. Use your job manager logs wisely to look for pain points, and modularize and parallelize when possible!
