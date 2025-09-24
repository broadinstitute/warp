# Changelog Style Guide Overview
The style guide provides formatting guidelines and language suggestions for pipeline changelogs. ...uggestions for pipeline changelogs. It is divided into four sections: 1) Changelog Information, which details the types of changes listed in a changelog, 2) Language Usage, which describes language style and syntax for changelog information, 3) Changelog Format, which provides sample formatting for changelog information, and 4) Sample Changelog Entries, which displays two examples of changelog entries taken from the [Optimus.changelog.md file](https://github.com/broadinstitute/warp/blob/develop/pipelines/wdl/optimus/Optimus.changelog.md).  

## Changelog Information
Pipeline changelogs have four informational parts:

#### 1) The changelog.md file name
The changelog file name should be consistent across pipelines. The format is "<pipeline name>.changelog.md".
 *   Ex: Optimus.changelog.md

The file name is not included in the actual changelog- it is just the file name.

#### 2) The pipeline release version name 
  *  Ex: [Illumina Genotyping Array pipeline](https://github.com/broadinstitute/warp/blob/develop/pipelines/wdl/genotyping/illumina/IlluminaGenotypingArray.changelog.md) version name is a number (i.e. "1.0.0")

The version name is listed at the top of each changelog entry section. It should follow [Semantic versioning (SemVer)](https://semver.org/) which uses the major.minor.patch version number. 

#### 3) The date of last commit (YYYY-MM-DD)

The purpose of the date of last commit is to inform users of the relative timing of changes to the pipeline. This is **not a release date**, which would have to be added to the changelog post-release. This date is updated as each change is created. For each pipeline release in the changelog, there will only be one final commit date (as opposed to a commit date for every change in that particular release). 
  
#### 4) Details of the version changes

This section describes (in bullet point format) every type of change made to the pipeline for the current pipeline version. Types of changes include: additions, updates, changes, deprecations, removals, fixes, and security updates. The type of change should be used as the first word of the bullet point (see examples below). These changes should be updated as pipeline changes are made even if the updated pipeline version has not been released. This will enable end-users to see what changes are coming when the new release is published.

If changes are breaking changes to the pipeline (a major version update), this changelog section will be divided into two subsections: "Breaking changes" and "Additional changes". An example of these headers is shown in the [Sample Changelog Entries](#sample-changelog-entries) section.

## Language Usage

All changelog notes should be bulleted (no period at the end of the final sentence of any one bullet point). Each bullet should include one type of change, but more than one sentence can be used to describe the individual change. Bullets should be written in complete sentences, but with the caveat that personal pronouns (“I”, “We”, “They”) are omitted. The first word of each bullet should be capitalized and describe the type of change made (in past tense). 

Examples of bullet points:

*  "Updated the WDL with the latest docker image”
*  "Added an output to the WDL. Users will now see this output after running the WDL"

For all bullet points, use active voice instead of passive voice. Examples are shown below. 

Examples of passive voice:
*  "Broken link in the StarAlign task **was fixed**"

Example of active voice:
*  "**Fixed** the broken link in the StarAlign task"

You can find more examples of active voice from the [University of Wisconsin-Madison's Writing Center](https://writing.wisc.edu/handbook/style/ccs_activevoice/).


## Changelog Format
The following is the markdown format required for all major changelog version updates which have breaking changes:

\# "insert Release Version Name here"

"insert date of last commit in YYYY-MM-DD here" (Date of Last Commit)

\### Breaking changes

* "insert description of first breaking pipeline change here" 
* "insert description of additional breaking changes here- add more bullets as necessary for additional changes"

\### Additional changes

* "insert description of first non-breaking pipeline change here" 
* "insert description of additional change here- add more bullets as necessary for additional changes"


The following is the markdown format required for non-major changelog version updates with non-breaking changes:

\# "insert Release Version Name here"

"insert date of last commit in YYYY-MM-DD here" (Date of Last Commit)

* "insert description of first pipeline change here" 
* "insert description of additional change here- add more bullets as necessary for additional changes"

## Sample Changelog Entries 

### Major Version Update (Breaking Change)

# 4.0.0

2020-08-10 (Date of Last Commit)

### Breaking changes
* Changed sample_id to input_id

### Additional changes 
* Added input_name as an optional input for user provided sample_id
* Passed pipeline_version to output loom file  
* Added input_id_metadata_field and input_name_metadata_field as optional input


### Non-major Version Update (Non-breaking Change)

# 1.4.0

2019-11-08 (Date of Last Commit)

* Added support for V3 chemistry
* Updated the documentation with additional information for the README and additional files for Loom schema and BAM tags
* Updated the Zarr output

## Syntax issues

 Not all valid markdown can be parsed and displayed in the release notes on GitHub. The following are known syntax issues:
 
 * code snippets, specified with the \` character are not supported (even when escaped)
 * double quotes need to be escaped: \\"
