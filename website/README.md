# WARP documentation Website

This documentation website is built using [Docusaurus 2](https://v2.docusaurus.io/), a modern static website generator.

### Installation

You need to install `yarn` before starting to work with this site. Run `brew install yarn` on macOS, and then:

```
$ yarn
```

### Local Development

```
$ yarn start
```

This command starts a local development server and open up a browser window. Most changes are reflected live without having to restart the server.

### Build

```
$ yarn build
```

This command generates static content into the `build` directory and can be served using any static contents hosting service.


### Deployment

```
$ GIT_USER=<Your GitHub username> USE_SSH=true yarn deploy
```

If you are using GitHub pages for hosting, this command is a convenient way to build the website and push to the `gh-pages` branch.

_We are working on a more automatic Github Actions based deployment mechanism_
