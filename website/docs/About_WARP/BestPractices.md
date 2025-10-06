---
sidebar_position: 1
---


# Best Practices for Building Data Processing Pipelines

WARP pipeline development is guided by the best practices detailed below. We describe each of these best practices to give insight as to why they are important and we provide examples to give you a sense of how to apply them.

Overall, the best pipelines should be:

- automated
- easily testable
- portable
- scalable to their data
- easy to maintain

## Automation

### What is Automation?

Automation refers to the ability of a pipeline to run, end-to-end, without human intervention.

### Why do we care about automation?

Pipelines cannot scale to large amounts of data, or many runs, if manual steps must be performed within the pipeline. They also cannot be part of an automated system if they in fact are not automated. Manual steps will bottleneck your entire system and can require unmanageable operations. Moreover, manual steps performed by humans will vary, and will promote the production of data that can not be appropriately harmonized.

- _Do_: Reduce parameterization to minimal inputs that do not vary for each input data.
- _Do_: Remove the need for parameters, replacing them with data-driven settings.
- _Do_: Offer defaults that are generally applicable for inputs that cannot be defined in a data-driven manner.
- _Do_: Offer the ability to check the status of pipeline runs.
- _Don’t_: Assume any file produced at any step of the pipeline is ok. Always check the status of underlying tools (Eg. check return codes).
- _Don’t_: Keep output files produced by steps of the pipeline that errored; people will accidently use them if they exist. (Do keep logs for debugging.)
- _Don’t_: Delete outputs from steps that passed when the full pipeline fails, keeping them enables you to pick up where you left off.
- _Don’t_: Use tools that are “buggy” or fragile, find alternatives or improve the tools.

## Testability

### What is a testable pipeline?

A testable pipeline is one in which isolated sections or the full pipeline can be checked for specified characteristics without modifying the pipeline’s code. Testability requires the existence of appropriate data with which to run the test and a testing checklist that reflects a clear understanding of how the data will be used to evaluate the pipeline.

### Why do we care about testabilty?

The availability of test data enables validation that the pipeline can produce the desired outcome. Formulation of a testing checklist allows the developer to clearly define the capabilities of the pipeline and the parameters of its use.

- _Do_: Provide example test data with your pipeline/tool.
- _Do_: Provide the results of an execution of your pipeline/tool on the test data.
- _Do_: Refer to at least one real data set appropriate for your tool/pipeline with example output from an execution of your pipeline or tool.
- _Do_: Provide a checker tool.
- _Do_: Include an automated testing suite for pipeline tasks.
- _Do_: Include automated tests for the pipeline as an integrated unit (pipeline benchmarking).
- _Don’t_: Use tools that do not have automated testing suites.
- _Don’t_: Write tests that assume a specific instance of data.

## Portability

### What is pipeline portability?

Pipeline portability refers to the ability of a pipeline to execute successfully on multiple technical architectures.

### Why do we care about portability?

_Science._ Science is not science if results are not reproducible; the scientific method cannot occur without a repeatable experiment that can be modified. Data processing pipelines are an essential part of some scientific inquiry and where they are leveraged they should be repeatable to validate and extend scientific discovery.

_Impact._ Pipelines will have greatest impact when they can be leveraged in multiple environments. The more technical requirements for installing and running of a pipeline, the longer it will take for a researcher to have a usable running pipeline.

_Maintainability._ Over the long term, it is easier to maintain pipelines that can be run in multiple environments. Portability avoids being tied to specific infrastructure and enables ease of deployment to development environments.

To ensure that others will be able to use your pipeline, avoid building in assumptions about environments and infrastructures in which it will run.

### Configurability for running on different technical infrastructures.

Code should not change to enable a pipeline to run on a different technical architecture; this change in execution environment should be configurable outside of the pipeline code.

- _Do_: Use a workflow language that allows a separation between the code that executes the data processing logic and the logic to run the pipeline on an infrastructure. [WDL](https://software.broadinstitute.org/wdl/documentation) and [CWL](https://www.commonwl.org/user_guide/rec-practices) languages have this feature.
- _Don’t_: Put logic to run the pipeline in the same code that executes the logic to process the data.

### Separation between the environment, the execution of the pipeline, and the pipeline itself.

- _Do_: Use a containerization technology, such as [Docker](https://www.docker.com/), to execute software.
- _Do_: Incorporate into your testing the execution of the pipeline in multiple execution environments.
- _Don’t_: Put environmental paths in software tools or the workflow language. When they must exist they belong in the configuration or (if they refer to the execution environment) in the container’s build instructions (Eg. Dockerfile).

## Scaling Characteristics

### What do we mean by scaling characteristics?

Scaling characteristics describe the performance of the pipeline given a certain amount of data. This is often described with Big O notation when describing algorithms. This answers the question: As the size of the data for the pipeline increases, how many additional computes are needed to process that data? One would want to avoid algorithms or tools that scale poorly, or improve this relationship to be linear (or better).

### Why care about scalability?

If you have poor scaling characteristics, it may take an exponential amount of time to process more data. This will eventually require unreasonable amounts of time (and money if running in the cloud) and generally reduce the applicability of the pipeline.

- _Do_: Measure the relationship between the size of your input (or something equivalent) and resources needed to successfully complete the pipeline.
- _Don’t_: Incorporate tools that have poor scaling characteristics unless they offer significant improvement to the pipeline.

## Maintainability

### What is a maintainable pipeline?

A pipeline that can be easily operated and updated is maintainable.

### Why do we care about maintainability?

The majority of the life of code involves maintenance and updates. Design and initial implementation require vastly shorter amounts of time compared to the typical time period over which the code is operated and updated. This is generally true in many areas of software engineering. Thus it is important to engineer software so that the maintenance phase is manageable and does not burden new software development or operations.

### Readability

Software is a living document that should be easily read and understood, regardless of who is the reader or author of the code.

- _Do_: Work in a space where the code is easy to access and navigate (Eg. [GitHub](https://github.com/))
- _Do_: Use common software package structure and idioms to aid the navigation of the software.
- _Do_: Use automated documentation for all technical documents as much as possible.
- _Don’t_: Write a large amount of documentation that does not live beside or within the code itself (it will become out of date).

### Modularity

Modularity enables small units of code to be independently benchmarked, updated, validated, and exchanged as science or technology changes. Using these small units enables more rapid updates and better adaptation to innovation.

- _Do_: Save progress by creating intermediate output between modules as they successfully complete.
- _Don’t_: Make monolithic tasks that perform many functionalities for the sake of speed.
- _Don’t_: Break every functionality of a pipeline into a separate module. (This contrasts with not making monolithic tasks; there is an optimum between monolithic tasks and highly resolved modularity that is the goal. One can use benchmarking, the tendency for functionality to be updated, and how dependent functionalities are to get a sense of what should be separate and what can be combined.)

### Leveraging Standards

We recommend using standard file formats and interfaces. In computational biology, [GA4GH](https://www.ga4gh.org/genomic-data-toolkit/) is a great source of these standards. In cases where new formats are needed, we recommend working with a standards group like [GA4GH](https://www.ga4gh.org/) if possible.

- _Do_: When using containerization technologies, follow best practices to assure associated images do not update without explicit updates.
- _Do_: Make both the images and the build files (Dockerfile) available to document the environment. [More on Dockerfiles](https://docs.docker.com/develop/develop-images/dockerfile_best-practices/).

## Cross-Platform Compatibility
This section offers resources to help you adapt WDLs for use across various environments. The links below point to key tools and platforms that support WDL execution and provide guidance for modifying workflows to fit specific contexts.


 **1. Cromwell Documentation**

Cromwell is a popular WDL execution engine that supports running workflows in various environments. Refer to the official [Cromwell documentation](https://cromwell.readthedocs.io/en/stable/) for guidance on adapting WDLs for use in local, cloud, or hybrid environments.

 **2. MiniWDL**
 
MiniWDL is a lightweight version of WDL designed for local or smaller-scale workflows. It can help in environments where a full Cromwell installation may not be necessary. Refer to the official [MiniWDL documentation](https://miniwdl.readthedocs.io/en/latest/) for more details.

 **3. Seven Bridges Genomics (Velsera)**
 
[Seven Bridges Genomics](https://www.sevenbridges.com/), now known as Velsera, provides a cloud platform for bioinformatics workflows. The platform supports WDL and provides tools for modifying workflows to run on their infrastructure.

 **4. Manifold**
 
[Manifold](https://www.manifold.ai/) is a cloud-native platform that supports WDL, CWL, and Nextflow, enabling scalable and reproducible bioinformatics workflows. It offers tools and documentation to help adapt WDLs for efficient execution across cloud and hybrid environments.

## Versioning

Versioning pipelines and associated Docker images allows you to determine when and how data is created (provenance). As you make improvements and changes to your pipeline, it is important to know which version of the pipeline and software you used to create a given dataset so that it can be easily reproduced. This not only facilitates scientific reproducibility for the greater community, it also allows you to verify that new pipeline changes produce consistent results. We recommend choosing a consistent versioning system (for example, the [semantic system](https://semver.org/)) and tracking pipeline changes in a [changelog](https://keepachangelog.com/en/1.0.0/).

## Licensing

### What is licensing?

According to Wikipedia "A software license is a legal instrument (usually by way of contract law, with or without printed material) governing the use or redistribution of software.” (see [this Wikipedia article](https://en.wikipedia.org/wiki/Software_license) for details).

### Why do we care about licensing?

:::warning NOTE
This section is opinion and is NOT legal advice.
:::

Licenses sometimes legally bind you as to how you use tools, and sometimes the terms of the license transfer to the software and data that is produced. This can restrict the potential for leveraging the pipeline and may require additional work.

- _Do_: Select tools that are openly licensed to run in your pipelines to avoid the possibility that legal requirements will restrict execution where technical requirements do not.
- _Don’t_: Create software tools or libraries without licenses, clear guidance on your intent for use is important.
