/** @type {import('@docusaurus/types').DocusaurusConfig} */
const path = require("path");

module.exports = {
  title: 'WARP',
  tagline: 'Warp Analysis Research Pipelines',
  url: 'https://broadinstitute.github.io',
  baseUrl: '/warp/', // FIXME: /warp/
  onBrokenLinks: 'warn',
  onBrokenMarkdownLinks: 'warn',
  favicon: 'img/favicon.ico',
  trailingSlash: false,
  organizationName: 'broadinstitute', // Usually your GitHub org/user name.
  projectName: 'warp', // Usually your repo name.
  plugins: [
    path.resolve(__dirname, 'plugins/docusaurus-plugin-hotjar')
  ],
  themeConfig: {
    docs: {
      sidebar: {
        hideable: true,
      },
    },
    //  sidebar.hideable: true,
    navbar: {
      title: 'WARP Pipelines',
      // logo: {
      //   alt: 'My Site Logo',
      //   src: 'img/logo.svg',
      // },
      items: [
        {
          type: 'doc',
          docId: 'get-started',
          position: 'right',
          label: 'Documentation',
        },
        {to: '/blog', label: 'Blog', position: 'right'},
        {
          href: 'https://github.com/broadinstitute/warp',
          position: 'right',
          className: 'header-github-link',
          'aria-label': 'GitHub repository',
        },
      ],
    },
    // to integrate Hotjar feedback
    hotjar: {
      siteId: '2427684',
    },
    footer: {
      style: 'dark',
      links: [
        {
          title: 'Docs',
          items: [
            {
              label: 'Get Started',
              to: '/docs/get-started',
            },
          ],
        },
        {
          title: 'Guide and Policy',
          items: [
            {
              label: 'Contribution Guide',
              href: '/docs/contribution/README',
            },
            {
              label: 'Privacy',
              href: '/privacy',
            },
          ],
        },
        {
          title: 'Resources',
          items: [
            {
              label: 'Blog',
              to: '/blog',
            },
            {
              label: 'Changelog',
              to: 'https://github.com/broadinstitute/warp/releases',
            },
            {
              label: 'GitHub',
              href: 'https://github.com/broadinstitute/warp',
            },
          ],
        },
      ],
      copyright: `Copyright © ${new Date().getFullYear()} Copyright © Data Sciences Platform, Broad Institute.`,
    },
    algolia: {
      appId: 'BX0CVVQ6KQ',
      indexName: 'broadinstitute',
      apiKey: 'a724a683820e05c85cbec97fcd24819d',
      contextualSearch: false,
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        gtag: {
          trackingID: 'G-NZ4BLRHQYS',
          anonymizeIP: true,
        },
        docs: {
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          editUrl:
            'https://github.com/broadinstitute/warp/edit/develop/website/',
          showLastUpdateAuthor: true,
          showLastUpdateTime: true,
          sidebarCollapsible: true,
        },
        blog: {
          showReadingTime: true,
          // Please change this to your repo.
          editUrl:
            'https://github.com/broadinstitute/warp/edit/develop/website/blog/',
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
};
