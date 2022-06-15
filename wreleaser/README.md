<p align="center">
        <img src="assets/wreleaser-logo.png">
</p>

<p align="center">
        <i><b>wreleaser</b> (WARP releaser)</i> is  a simple CLI tool for querying WARP <a href="https://github.com/broadinstitute/warp/releases">releases<a>.
</p>

## :dna: Overview

WARP maintains and releases over two dozen cloud-optimized pipelines for biological data processing. Finding information for various pipeline releases via the repo's web interface is clunky and tedious.

_wreleaser_ serves as a wrapper around the [GitHub Releases API](https://docs.github.com/en/rest/reference/repos#releases) specifically for WARP releases.

### How does it work?

On initial invocation _wreleaser_ will hit the GitHub Releases API for WARP and cache all of the releases in `$HOME/.wreleaser/cache.json`.

WARP doesn't release frequently so to avoid having to hit the API regularly _wreleaser_ will return results from the cache by default.

To reset your cache you can run `wreleaser reset`.
## :gear: Installation

#### Option 1

Make sure you have Go installed ([download](https://jimkang.medium.com/install-go-on-mac-with-homebrew-5fa421fc55f5)). Then install _wreleaser_ with the `go get` command:

```bash
go install github.com/broadinstitute/warp/wreleaser@latest
```

#### Option 2

We have binaries available [here](./scripts) for `linux/amd64` and `darwin/amd64`.

If you have a different OS/Architecture please clone this repo and run [`go build`](https://pkg.go.dev/go/build) from this directory to compile the binary for your system.

Once the binary is on your machine please add it to a location in your `$PATH`.
## :star: Quickstart

```bash
$ wreleaser --help

$ wreleaser is a lightweight CLI to list, query and download various WARP releases

Currently available pipelines:
        - AnnotationFiltration
        - Arrays
        - BroadInternalRNAWithUMIs
        - BroadInternalArrays
        - BroadInternalImputation
        - BroadInternalUltimaGenomics
        - CEMBA
        - CheckFingerprint
        - CramToUnmappedBams
        - ExomeGermlineSingleSample
        - ExomeReprocessing
        - ExternalExomeReprocessing
        - ExternalWholeGenomeReprocessing
        - GDCWholeGenomeSomaticSingleSample
        - IlluminaGenotypingArray
        - Imputation
        - JointGenotyping
        - JointGenotypingByChromosomePartOne
        - JointGenotypingByChromosomePartTwo
        - MultiSampleArrays
        - MultiSampleSmartSeq2
        - MultiSampleSmartSeq2SingleNucleus
        - Optimus
        - ReblockGVCF
        - RNAWithUMIsPipeline
        - SmartSeq2SingleNucleus
        - SmartSeq2SingleSample
        - UltimaGenomicsWholeGenomeGermline
        - UltimaGenomicsJointGenotyping
        - ValidateChip
        - VariantCalling
        - WholeGenomeGermlineSingleSample
        - WholeGenomeReprocessing

Usage:
  wreleaser [command]
```

## :eyes: Examples

### Commands

#### `wreleaser info all`
Displays release information for all pipelines

Flags:

* -l, --latest  Retrieve only the latest releases for each pipeline
* -h, --help    Display help text
* --output      Pipe stdout to file

```bash
$ wreleaser info all # Gets all the releases for every pipeline

$ wreleaser info all --latest # Gets only the latest releases for every pipeline

$ wreleaser info all --latest --output=foo.json # Gets the latest releases for every pipeline and write it to foo.json
```

#### `wreleaser info list`
Display release information for one or more pipelines

Flags:

* -l, --latest  Retrieve only the latest releases for each pipeline
* -v, --version Specific version number of pipeline to retrieve (can only use for single pipeline arguments)
* -h, --help    Display help text
* --output      Pipe stdout to file

```bash
$ wreleaser info list Arrays # Gets all the releases for the Arrays pipeline

$ wreleaser info list Arrays --version=v2.3.1 # Gets the release for Arrays pipeline version v2.3.1

$ wreleaser info list Arrays ExomeGermlineSingleSample # Gets all releases for Arrays and ExomeGermlineSingleSample pipelines

$ wreleaser info list CEMBA Optimus --latest # Gets the latest releases for CEMBA and Optimus pipelines

$ wreleaser info list CEMBA Optimus --latest --output=foo.json # Gets the latest releases for CEMBA and Optimus pipelines and write it to foo.json
```

#### `wreleaser reset`
Clear the directory that caches WARP releases

Default directory is `$HOME/.wreleaser/cache.json` but you can choose to write/clear cache with the `--cachedir` flag for any command

Flags:

* h, --help Display help text
* --cachedir Directory to cache release information

```bash
$ wreleaser reset # Clears the cache at $HOME/.wreleaser/cache.json

$ wreleaser reset --cachedir=/path/to/custom # Clears the cache at /path/to/custom
```

## :boom: Troubleshooting
If you have any questions about this tool or need some help configuring it please reach out at [dsde-engineer@broadinstitute.org](dsde-engineer@broadinstitute.org).



