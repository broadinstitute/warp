# (WIP) WARP Website

A static website built with [VuePress](https://vuepress.vuejs.org/).

## Prerequisites

You need to install the following before setting up this project:

- `yarn`

On macOS, you can install the above by:

```bash
brew install yarn
```

## Setup

Install all of the dependencies by:

```bash
yarn install
```

## Local Development

In order to spin up the dev server locally for development:

```bash
yarn dev
```

### Routing

The routing table of the entire site can be found at [src/.vuepress/config.js](./src/.vuepress/config.js). See [VuePress Docs](https://vuepress.vuejs.org/guide/directory-structure.html#directory-structure) for more details.

### Custom Styles

Most of the customization is done in 3 places:

- The plugin section of `src/.vuepress/config.js`.
- `src/.vuepress/style/index.styl`.
- `src/.vuepress/style/palette.styl`.

### Build

To build the static site, from the root, run:

```bash
yarn build
```

this will generate a `dist` folder in `src/.vuepress`.

### Formatting

To format with default code styles, run:

```bash
yarn format
```
(this will only format files under `.vuepress`, you could change the `"format": "prettier --write src"` in `[package.json](./package.json)` to also format markdown file if needed)
