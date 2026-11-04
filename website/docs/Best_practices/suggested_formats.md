---
sidebar_position: 5
---

# WDL formatting tips 
There are multiple WDL resources on how to keep scripts clean and easy to read. The Broad Pipeline Development team suggests checking out the [BioWDL style guide](https://biowdl.github.io/styleGuidelines.html) for tips and tricks. The team follows most of the BioWDL guidelines, with a few exceptions and additions detailed in the table below:

| Style attribute | BioWDL recommendation | Pipeline Development team WDL practice | 
| --- | --- | --- |
| Indentation spacing | [Indent 4 spaces](https://biowdl.github.io/styleGuidelines.html#1-indentation) | Indent 2 spaces. | 
| Blank lines | [Add blank lines between workflow sections](https://biowdl.github.io/styleGuidelines.html#2-blank-lines). | No particular convention for blank lines. In general, do not add blank lines between the closing braces of a parent and child block. |
| Line breaks | [Adhere to line break limits](https://biowdl.github.io/styleGuidelines.html#4-line-length-and-line-breaks). | Do not follow specific formatting for line breaks. |
| Naming conventions | [Use the listed naming standards](https://biowdl.github.io/styleGuidelines.html#5-naming-conventions). | Do not follow specific naming conventions. In general, variables are lowercase underscore (python style) and alias calls use UpperCamelCase. |
| Task command section input arguments | N/A. No guidelines provided on formatting input arguments. | Use one input argument per code line. For an example, see the STAR parameters in the [StarAlign task WDL](https://github.com/broadinstitute/warp/blob/develop/tasks/wdl/StarAlign.wdl). |
| Filenames | N/A. No guidelines provided. | Recommends renaming files to meet tool requirements, including adding or removing necessary file extensions. For example, some tools require a “.gz” extension). |

## Resources
* [BioWDL Styleguide](https://biowdl.github.io/styleGuidelines.html)
