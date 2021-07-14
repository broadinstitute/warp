# wreleaser

**_wreleaser_** (_WARP releaser_) is a simple CLI tool for querying WARP [releases](https://github.com/broadinstitute/warp/releases)

## :dna: Overview

WARP maintains and releases over a dozen cloud-optimized pipelines for biological data processing. Finding information for various pipeline releases via the repo's web interface is clunky and tedious.

_wreleaser_ serves as a wrapper around the [GitHub Release API](https://docs.github.com/en/rest/reference/repos#releases) specifically for WARP releases.
## :gear: Installation

#### Option 1

Make sure you have Go installed ([download](https://jimkang.medium.com/install-go-on-mac-with-homebrew-5fa421fc55f5)). Then install _wreleaser_ with the `go get` command:

```bash
go get -u github.com/broadinstitute/warp/wreleaser@wd_GL-1606-warp-releases
```

#### Option 2

We have binaries available [here](./scripts) for `linux/amd64` and `darwin/amd64`.

If you have a different OS/Architecture please clone this repo and run [`go build`](https://pkg.go.dev/go/build) to compile the binary for your system.

Once the binary is on your machine please add it your `$PATH`.
## :star: Quickstart

```bash
$ wreleaser --help
```
