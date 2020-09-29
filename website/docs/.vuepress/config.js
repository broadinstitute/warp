const { description } = require("../package");

module.exports = {
  base: "/warp/",
  /**
   * Ref：https://v1.vuepress.vuejs.org/config/#title
   */
  title: "WARP",
  /**
   * Ref：https://v1.vuepress.vuejs.org/config/#description
   */
  description: description,

  /**
   * Extra tags to be injected to the page HTML `<head>`
   *
   * ref：https://v1.vuepress.vuejs.org/config/#head
   */
  head: [
    ["meta", { name: "theme-color", content: "#3eaf7c" }],
    ["meta", { name: "apple-mobile-web-app-capable", content: "yes" }],
    [
      "meta",
      { name: "apple-mobile-web-app-status-bar-style", content: "black" },
    ],
    ["link", { rel: "icon", href: "/logo.png" }],
    ["meta", { name: "og:locale", content: "en_US" }],
  ],

  /**
   * Theme configuration, here is the default theme configuration for VuePress.
   *
   * ref：https://v1.vuepress.vuejs.org/theme/default-theme-config.html
   */
  themeConfig: {
    repo: "https://github.com/broadinstitute/warp",
    editLinks: false,
    docsDir: "",
    sidebarDepth: 3,
    editLinkText: "",
    smoothScroll: true,
    editLinks: true,
    lastUpdated: true,
    nav: [
      {
        text: "Documentation",
        link: "/documentation/",
      },
      {
        text: "Blog",
        link: "/blog/",
      },
      {
        text: "Learn More",
        items: [
          {
            text: "Explore",
            items: [{ text: "Research", link: "/research/" }],
          },
          {
            text: "Resources",
            items: [
              { text: "Contribution Guide", link: "/contribution/" },
              {
                text: "Changelog",
                link: "https://github.com/broadinstitute/warp/releases",
              },
            ],
          },
        ],
      },
    ],
    sidebar: {
      "/documentation/": [
        {
          title: "BestPractices",
          collapsable: true,
          children: ["BestPractices"],
        },
        {
          title: "Optimus",
          collapsable: false,
          children: ["optimus/", "optimus/optimus.methods"],
        },
      ],
    },
  },

  /**
   * Apply plugins，ref：https://v1.vuepress.vuejs.org/zh/plugin/
   */
  plugins: [
    // "@vuepress/active-header-links",
    "@vuepress/plugin-back-to-top",
    "@vuepress/plugin-medium-zoom",
    ["vuepress-plugin-mathjax", { target: "svg", macros: { "*": "\\times" } }],
  ],
};
