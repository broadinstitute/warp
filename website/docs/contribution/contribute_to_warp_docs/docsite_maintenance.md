---
sidebar_position: 1
---

# WARP Documentation Maintenance Guide

This documentation site is built on top of the `React`-based  framework [Docusaurus](https://docusaurus.io), so most of the recommended configuration can be found in the frameworks documentation.

Due to the nature of the frontend framework, most of the plugins this docsite is using can be deprecated or updated over the time. Github has the nice `dependabot` feature which periodically scans this code base and finds the potential vulnerbilities. There are usually 2 types of updates:

- Major version nump
- Security patches

The former one requires more human engagement and you should follow the official [Docusaurus](https://docusaurus.io) instructions about how to upgrade the major version. This will happen when docusaurus is in GA phase, and perhaps when it has a major/minor version bump.

The latter one in most cases will be identified by `dependabot` and only requires a human review + approval, before commentting `@dependabot merge` on the PR. In rare cases, the PR opened by `dependabot` will have merge conflicts against the `develop` branch. Developers will need to `cherry-pick` the commit made by `dependabot`, resolve the conflicts and open a separate PR for it.

When you want to test your changes to the docs, you can deploy the site locally.
If doing this for the first time, install the dev version of the site using the root of the repo:
```
yarn --cwd=website install
````
Then deploy the dev version using:

```
yarn --cwd=website start
```

to spin up a local server to preview your changes.

For more details on deployment, check this [page](https://docusaurus.io/docs/deployment).
