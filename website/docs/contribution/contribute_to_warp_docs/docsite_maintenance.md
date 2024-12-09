---
sidebar_position: 1
---

# WARP Documentation Maintenance Guide

This documentation site is built on top of the `React`-based  framework [Docusaurus](https://docusaurus.io), so most of the recommended configuration can be found in the frameworks documentation.

## Versions and dependabot

Due to the nature of the frontend framework, most of the plugins this docsite is using can be deprecated or updated over the time. Github has the nice `dependabot` feature which periodically scans this code base and finds the potential vulnerbilities. There are usually 2 types of updates:

- Major version nump
- Security patches

The former one requires more human engagement and you should follow the official [Docusaurus](https://docusaurus.io) instructions about how to upgrade the major version. This will happen when docusaurus is in GA phase, and perhaps when it has a major/minor version bump.

The latter one in most cases will be identified by `dependabot` and only requires a human review + approval, before commentting `@dependabot merge` on the PR. In rare cases, the PR opened by `dependabot` will have merge conflicts against the `develop` branch. Developers will need to `cherry-pick` the commit made by `dependabot`, resolve the conflicts and open a separate PR for it.

## Deploying local site for testing and editing

When you want to test your changes to the docs, you can deploy the site locally.
* first you need yarn, which is installed with npm: [mac install](https://classic.yarnpkg.com/lang/en/docs/install/#mac-stable)
* then you need to build the local website in the root of the docs directory `warp/website` on your desired branch
   * `yarn --cwd=website install` # install the site
   * `yarn --cwd=website start` # spin up a local server to preview your changes
   * the local website should popup automatically in your default browser
   * the website should update automatically as you edit the markdown files

For more details on deployment, check this [page](https://docusaurus.io/docs/deployment).
