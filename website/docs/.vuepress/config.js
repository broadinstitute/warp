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
        link: "/documentation/Intro",
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
          title: "About_WARP",
          collapsable: false,
          children: [
            "About_WARP/BestPractices",
            "About_WARP/PipelineRequirements",
            "About_WARP/TestingPipelines",
            "About_WARP/VersionAndReleasePipelines"
          ],
        },
        {
          title: "Pipelines",
          collapsable: false,
          children: [
            {
              title: "CEMBA MethylC Seq Pipeline",
              collapsable: true,
              children: [
                "Pipelines/CEMBA MethylC Seq Pipeline/CEMBA.methods",
                "Pipelines/CEMBA MethylC Seq Pipeline/"
              ],
            },
            {
              title: "Single Cell ATAC Seq Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Single Cell ATAC Seq Pipeline/",
                "Pipelines/Single Cell ATAC Seq Pipeline/scatac.methods"
              ],
             },
            {
              title: "Exome Germline Single Sample Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Exome Germline Single Sample Pipeline/"
              ],
            },
            {
              title: "Illumina Genotyping Arrays Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Illumina Genotyping Arrays Pipeline/IlluminaGenotypingArray.documentation"
              ],
            },
            {
              title: "Optimus Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Optimus Pipeline/",
                "Pipelines/Optimus Pipeline/Bam_tags",
                "Pipelines/Optimus Pipeline/Loom_schema",
                "Pipelines/Optimus Pipeline/optimus.methods"
              ],
            },
            {
              title: "Smart-seq2 Multi Sample Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Smart-seq2 Multi Sample Pipeline/",
                "Pipelines/Smart-seq2 Multi Sample Pipeline/Loom_schema",
                "Pipelines/Smart-seq2 Multi Sample Pipeline/smart-seq2.methods"
              ],
            },
            {
              title: "Smart-seq2 Single Sample Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Smart-seq2 Single Sample Pipeline/",
                "Pipelines/Smart-seq2 Single Sample Pipeline/Loom_schema"
              ],
            },
            {
              title: "Whole Genome Germline Single Sample Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Whole Genome Germline Single Sample Pipeline/"
              ],
            },
          ],
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
