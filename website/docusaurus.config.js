/** @type {import('@docusaurus/types').DocusaurusConfig} */
const path = require("path");

module.exports = {
  title: 'WARP',
  tagline: 'WDL Analysis Research Pipelines',
  url: 'https://broadinstitute.github.io',
  baseUrl: '/warp/', // FIXME: /warp/
  onBrokenLinks: 'warn',
  onBrokenMarkdownLinks: 'warn',
  favicon: 'img/favicon.ico',
  organizationName: 'broadinstitute', // Usually your GitHub org/user name.
  projectName: 'warp', // Usually your repo name.
  plugins: [
    path.resolve(__dirname, 'plugins/docusaurus-plugin-hotjar')
  ],
  themeConfig: {
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
    googleAnalytics: {
      trackingID: 'UA-181048107-1',
      // the Google Analytics Tracking ID
      anonymizeIP: true, // Should IPs be anonymized?
    },
    algolia: {
      apiKey: '57c6f6f64593f4fbec4c87a98d80f7e8',
      indexName: 'broadinstitute',
      contextualSearch: false,
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          editUrl:
            'https://github.com/broadinstitute/warp/edit/develop/website/',
          showLastUpdateAuthor: true,
          showLastUpdateTime: true,
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
