---
title: All of Us Pipelines
className: aou-doc-page
---

<div className="aou-folder-text">

<div className="aou-hero">
	<p>
		<strong>All of Us Analysis Workflows</strong><br />
		Reproducible pipelines for community review, reuse, and adaptation.
	</p>
</div>

This section contains workflows used in support of genomic analyses for the [All of Us Research Program](https://allofus.nih.gov/). We release these pipelines to support reproducibility and to make methods easier for the community to review, reuse, and adapt.

## Pipeline Overview

The following All of Us analysis areas are currently available or in-progress in the WARP repository:

* [Admixture](https://github.com/broadinstitute/warp/blob/develop/all_of_us/admixture/README.md)
* [Ancestry](https://github.com/broadinstitute/warp/blob/develop/all_of_us/ancestry/README.md)
* [CMRG (Challenging Medically Relevant Genes)](https://github.com/broadinstitute/warp/blob/develop/all_of_us/cmrg/README.md)
* [HLA Genotyping](https://github.com/broadinstitute/warp/blob/develop/all_of_us/hla/README.md)
* [Mitochondria](https://github.com/broadinstitute/warp/blob/develop/all_of_us/mitochondria/README.md)
* [PCA](https://github.com/broadinstitute/warp/blob/develop/all_of_us/PCA/README.md)
* [Phasing](https://github.com/broadinstitute/warp/blob/develop/all_of_us/phasing/README.md)
* [RNA-seq](https://github.com/broadinstitute/warp/blob/develop/all_of_us/rna_seq/README.md)

## Versioning and Changelogs

<div className="aou-callout">
	<strong>Version note:</strong> All of Us workflow versions (for example, <code>aou_9.0.0</code>) align to the release in which a pipeline was used.
</div>

Pipeline WDLs in this area include companion changelog files (for example, `*.changelog.md`). The All of Us versioning convention (for example, `aou_9.0.0`) reflects the release version in which the pipeline was used.

## Docker Images

<div className="aou-callout">
	<strong>Container note:</strong> Docker images for these workflows are maintained in the sister repository <a href="https://github.com/broadinstitute/warp-tools"><code>warp-tools</code></a>.
</div>

Several Docker images used by these workflows are maintained in the sister repository [`warp-tools`](https://github.com/broadinstitute/warp-tools). These Dockers were created to ensure consistency across pipeline runs.

## What to Expect

<div className="aou-callout">
	<strong>Community expectations:</strong> Some analyses are still in progress, pipelines are provided as-is, and suggestions are welcome.
</div>

* Some analyses in this folder are still in progress.
* Pipelines are provided as-is. These pipelines are designed to run in a Cromwell environment such as Terra.
* We welcome suggestions and improvement ideas from the community. Please post to warp-issues.

</div>
