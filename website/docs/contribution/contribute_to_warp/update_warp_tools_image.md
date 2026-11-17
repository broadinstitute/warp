---
sidebar_position: 3
---

# Updating a warp-tools Image

Most WARP pipelines run Docker images built in [warp-tools](https://github.com/broadinstitute/warp-tools). WARP pins each one as `<image>:<release-tag>@sha256:<digest>`: the digest fixes what runs, and the release tag says which release it is. A change takes two pull requests (PRs), one in warp-tools and one in WARP, and four steps.

## Workflow

The example is the scvi-scanvi 1.7.2 update for scANVI.

1. **Change the image in warp-tools and open a PR to `develop`.** Bump `DOCKER_IMAGE_VERSION`, and the tool version if it changed, in the image's `docker_build.sh`. CI builds the image and pushes it under the branch name.
2. **From your WARP branch, release and pin the build.** Once the PR's CI build has finished, pass the script the warp-tools PR (number, URL or branch name):

   ```bash
   python3 scripts/image_pins.py release 222
   # Tagged us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:… as 1.7.2-1.5.1-<timestamp>
   # pipelines/wdl/scanvi/scANVI.wdl:77: … -> …/scvi-scanvi:1.7.2-1.5.1-<timestamp>@sha256:…
   ```

   The script finds what to do from the PR: each image built by its latest commit, the build's digest, and the version from `docker_build.sh`. It refuses if the build hasn't finished or `DOCKER_IMAGE_VERSION` wasn't bumped. It then tags the build `<version>-<unix-timestamp>` and pins every WDL reference to the image; give paths after the PR to pin only some. Tagging needs write access to the registry (`gcloud auth`).
3. **Bump, test and review in WARP.** Bump each affected pipeline and its changelog (see the rules below), validate with womtool and `scripts/validate_release.sh -g origin/staging`, and run the tests: Plumbing, plus Scientific for a minor or major bump. If the image has to change again, push to the warp-tools PR and repeat step 2 once it has built; pushes that don't change the image need nothing.
4. **Merge the warp-tools PR, then the WARP PR.**

## Rules

- **Pins** are always `<image>:<release-tag>@sha256:<digest>`, written by the script. Never pin `latest`, a branch tag, a bare tag or a bare digest. The advisory `Image pins` check reports any that aren't.
- **Release tags never move.** A fix is a new release with a new tag. The script never retags; it reuses a build's existing tag for the same version.
- **Version bump:** minor by default, because a rebuilt image may shift outputs slightly; major when outputs change qualitatively; patch, with a "No functional impact" entry, for a pipeline that only imports an edited WDL and runs no task using the new image. See [Version and Release Pipelines](../../About_WARP/VersionAndReleasePipelines.md).
- **Changelog:** name the image, the release tag, what changed in the image and the effect on outputs, for example "Updated the scvi-scanvi image to 1.7.2-1.5.1-&lt;timestamp&gt;, which fixes how an AIT reference's counts are restored; outputs may differ slightly, with no reprocessing needed". Leave the digest to the WDL.
- **Shared verification pins:** the script skips `verification/VerifyTasks.wdl`, which is in about 27 test path filters; don't move it in a feature PR.
- **Docs pages** link to the WDL instead of copying image tags or digests.

## Outside the standard workflow

These are the only sanctioned ways to tag or pin differently. Release tags still never move, and a merged pin is still `<image>:<release-tag>@sha256:<digest>`.

| Case | What to do |
| --- | --- |
| No registry write access | Ask a maintainer who has it to run step 2. |
| Release an earlier build than the PR's latest, such as the one you already tested | Add `--digest sha256:<digest>` to step 2. |
| `docker_build.sh` computes part of the tag at build time (glimpse2, sv-imputation-rust-tools, zcall, illumina-iaap-autocall) | Add `--version <image-version>-<tool-version>` to step 2. |
| Pin an existing release | `python3 scripts/image_pins.py pin <image> <release-tag> [paths]` pins without tagging. |
| An image built outside a warp-tools PR | `gcloud container images add-tag us.gcr.io/broad-gotc-prod/<image>@sha256:<digest> us.gcr.io/broad-gotc-prod/<image>:<version>-$(date +%s)`, then `pin` that tag. |
| Test build under a custom tag | Run the image's warp-tools build workflow by hand with `image_tag`. Never name an existing tag; it would be overwritten. |
| Local builds | Use `docker build` for local testing. Never push a release-style tag from a laptop. |
| Azure | `crane copy us.gcr.io/broad-gotc-prod/<image>@sha256:<digest> dsppipelinedev.azurecr.io/<image>:<release-tag>` keeps the digest. |
