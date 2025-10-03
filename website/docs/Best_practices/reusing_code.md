---
sidebar_position: 4
---
# Reusing WDL code
Whether you’re setting up workflows for future reusability or trying to reuse someone else’s code, the Broad Pipeline Development team has a few tips and tricks to consider for your use case.

## Setting up your code for reusability
Set up your WDL code to reuse tools and software keeping in mind the following tips:

### Modularize reusable code into small WDL tasks
If you’re going to reuse software tools and commands across workflows, it makes sense to thoughtfully modularize these commands into separate WDL tasks that can be called independently. When thinking about how or what to modularize, keep maintenance in mind; when you share a WDL task across workflows, updating that task impacts any other workflow that uses it, including any automated scripts that validate your workflows. You want to keep tasks small and functional enough that you can reuse them, but not so small that you increase your maintenance by having to update multiple tasks with a single workflow update. 

### Code WDL tasks in external WDL scripts, not in workflow definitions
When coding extensive WDL workflows that reuse multiple tasks, it makes sense to place those tasks in separate WDL scripts that your workflow(s) can import. This keeps your workflows clean and easier to read. It also allows you to choose how you group tasks together for reusability and maintenance. For example, you might choose to group together tasks that are for a specific workflow or group tasks by function. In either case, you can put tasks into a separate WDL script dedicated to the tasks alone. WARP has several examples of these task WDLs that you can use as a template, such as the [Imputation task WDL](https://github.com/broadinstitute/warp/blob/master/tasks/wdl/ImputationTasks.wdl), which groups tasks that are for the Imputation pipeline, or the [Utilities task WDL](https://github.com/broadinstitute/warp/blob/master/tasks/wdl/Utilities.wdl), which contains common genomic sequencing file conversion and validation tasks that are applied across multiple workflows. 

## Reusing someone else’s WDL code
The best way to reuse WDL code written by someone else depends on how you want to implement and update it for your workflows. You can reuse code exactly as written by the contributor or you can modify it for your own purposes. If you plan to use code as-is, it can be helpful to pull that code directly from whatever repository stores it. For example, if it’s in GitHub, you can either fork the repository or point directly to code using raw, full-length `http` paths (see Cromwell’s documentation on imports). You can also use forking if you’re planning to modify code, but in that case, sometimes it might even be easier to copy and paste the code. No matter your intentions, keep the following tips in mind when using external code.

### Decide if you want to keep the code up-to-date with the source
Sometimes you might want to reuse code that is regularly updated by external contributors. If you want to keep up with these updates, you don’t want to deal with the hassle of manually making modifications. This is why the Pipeline Development team recommends importing the WDL code directly from whatever respository it lives in. If you’re only trying to reuse a specific version of WDL code, you can still pull it directly, but you’ll want to pin to the specific code version.

### When modifying code, ensure Docker images are compatible with code adjustments
Another use case for reusing code is simply avoiding reinventing the wheel. The WDL community is wide and there are multiple existing code snippets and tasks that have just the tool you need for your analysis. However, keep in mind that some of these tasks are set up to use parameters and scripts that are specific to the Docker specified in the task runtime. If you’re planning on modifying existing code but still using the tool in the task’s Docker, it’s important to check that your modifications are compatible. 
