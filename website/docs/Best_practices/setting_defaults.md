---
sidebar_position: 3
---

# Setting default values
In WDL, default values are those that allow your workflow to run in the absence of a user-defined attribute. These include primary inputs that you’ve hardcoded to a value, but they can also include inputs that are assigned to dynamically calculated functions ([like autosizing functions for setting runtime parameters](./autosize.md)). If you’re wondering whether a parameter has a default value or not, just ask yourself, “Will the workflow run if I don’t specify an attribute in the WDL’s input JSON?”

## Why set up defaults?
Default values can be used for a variety of reasons. Perhaps you want to hardcode a value because you don’t want others changing parameters while running the workflow. Or maybe you have standard values, like reference files, that you apply consistently across workflow runs regardless of the data you analyze. Overall, defaults can help ensure a workflow won’t keel over across different runs and user hands. 

Below, the Broad Pipeline Development team provides some best practices for building flexible workflows with default values.

## Tip 1: Decide what kinds of inputs you’re using
Workflows have multiple input types: samples or files that you want to analyze, parameters you want to use to set up runtime environments, dockers that contain your software, or inputs that are specific to the type of samples you're analyzing, like experimental parameters, information about species, age, etc. When thinking about setting defaults, you might want to consider what inputs change regularly (like sample files that you’re trying to analyze), or those that stay consistent (like reference data if you’re primarily using a particular species). In general, consistent values are good candidates to set to a default value.

## Tip 2: Decide if default values should be modifiable
While it makes sense to hardcode inputs to a default value if the inputs are consistent, you might still want to build in flexibility to allow those inputs to be modified. Consider scenarios where other workflow users might need to change a default value, but that change won’t necessarily break the workflow. For example, maybe you primarily run your workflow on mouse data using mouse reference files, but your workflow tasks and tools work equally well with human data. In this case, it makes sense to hardcode some input values to reflect mouse analysis, but also to provide flexibility for someone who wants to do human analysis.

In contrast, there are some inputs that might be necessary for the workflow to appropriately run, like runtime parameters that set up computational environments. For example, docker image inputs are one type of parameter that could break the pipeline if changed. You’ll often see docker file paths listed as configurable input in the workflow definition. If the wrong docker image is used, it might make the workflow fail, or worse, cause a security concern. These kinds of inputs should not be easily modified from an accompanying JSON file. When possible, the Pipeline Development team recommends that you hardcode paths to docker images into the task runtime attributes. 


## Tip 3: Know which inputs can be modified with a JSON
Now that we’ve considered our different input types, let’s review which inputs can be modified, meaning they can be specified with an accompanying input JSON. Modifiable inputs include those that are in the workflow definition’s input section AND inputs in a task’s input section when you allow for [nested inputs](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#computing-call-inputs). You can allow nested inputs by using the boolean flag `allowNestedInputs` in your workflow(s) definition’s meta section. This allows you to set up task-level inputs from your WDL’s input JSON.  While this feature works with Cromwell and in the bioinformatics platform Terra, you should note that it’s not fully supported, so use it with caution.

Inputs that can’t be modified are those that are outside of the workflow definition inputs section or those that are in the task inputs when nested inputs are not allowed. Keep this information in mind when you have inputs that you don’t want others to easily modify.

## Tip 4: Define modifiable defaults at the highest workflow level
Once you’ve considered your inputs and decided which ones should have default values, you’ll need to decide where in the workflow (the workflow definition or task) you want to specify the defaults. The Pipeline Development team recommends that you define defaults at the highest workflow level that makes sense, keeping in mind whether you want task inputs to be modifiable. For example, if you want to set a modifiable default value for an input reference file, assign the file to its file path in the workflow definition and pass the corresponding variable into the tasks using the workflow definition call sections. Using defaults at the highest level means that if you need to change the default value, you won’t have to do it for each individual task.

## Resources
* [Nested inputs description](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#computing-call-inputs)
