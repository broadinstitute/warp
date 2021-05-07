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

  markdown: {
    lineNumbers: false,
    extendMarkdown: (md) => {
      md.use(require("markdown-it-imsize"));
      md.use(require("markdown-it-footnote"));
      md.use(require("markdown-it-mark"));
    },
  },

  /**
   * Theme configuration, here is the default theme configuration for VuePress.
   *
   * ref：https://v1.vuepress.vuejs.org/theme/default-theme-config.html
   */
  themeConfig: {
    repo: "https://github.com/broadinstitute/warp",
    docsDir: "website/docs",
    sidebarDepth: 3,
    editLinkText: "Edit this page on GitHub",
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
              { text: "Privacy", link: "/privacy/" },
            ],
          },
        ],
      },
    ],
    sidebar: {
      "/documentation/": [
        {
          title: "Welcome to WARP",
          collapsable: true,
          children: [
            "Intro"
          ],
        },
        {
          title: "About WARP",
          collapsable: true,
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
                "Pipelines/CEMBA_MethylC_Seq_Pipeline/",
                "Pipelines/CEMBA_MethylC_Seq_Pipeline/CEMBA.methods"
              ],
            },
            {
              title: "Exome Germline Single Sample Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Exome_Germline_Single_Sample_Pipeline/"
              ],
            },
             {
              title: "Genomic Data Commons Whole Genome Somatic Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Genomic_Data_Commons_Whole_Genome_Somatic/"
              ],
            },           
            {
              title: "Illumina Genotyping Arrays Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Illumina_Genotyping_Arrays_Pipeline/IlluminaGenotypingArray.documentation",
                "Pipelines/Illumina_Genotyping_Arrays_Pipeline/Illumina_genotyping_array_spec"
              ],
            },
            {
              title: "Optimus Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Optimus_Pipeline/",
                "Pipelines/Optimus_Pipeline/optimus.methods",
                "Pipelines/Optimus_Pipeline/Bam_tags",
                "Pipelines/Optimus_Pipeline/Loom_schema"
              ],
            },
            {
              title: "Single Cell ATAC Seq Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Single_Cell_ATAC_Seq_Pipeline/",
                "Pipelines/Single_Cell_ATAC_Seq_Pipeline/scatac.methods"
              ],
             },
            {
              title: "Smart-seq2 Multi Sample Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Smart-seq2_Multi_Sample_Pipeline/",
                "Pipelines/Smart-seq2_Multi_Sample_Pipeline/smart-seq2.methods",
                "Pipelines/Smart-seq2_Multi_Sample_Pipeline/Loom_schema"
              ],
            },
            {
              title: "Smart-seq2 Single Sample Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Smart-seq2_Single_Sample_Pipeline/",
                "Pipelines/Smart-seq2_Single_Sample_Pipeline/Loom_schema"
              ],
            },
            {
              title: "Whole Genome Germline Single Sample Pipeline",
              collapsable: true,
              children: [
                "Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/",
                "Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/wgs.methods"                
              ],
            },
          ],
        },
      ],
      "/contribution/": [
        {
          title: "Contribution Guide",
          collapsable: false,
          children: [
            "",
            "changelog_style",
            "doc_style"
          ],
        },
      ]
    },
  },

  /**
   * Apply plugins，ref：https://v1.vuepress.vuejs.org/zh/plugin/
   */
  plugins: [
    "@vuepress/active-header-links",
    "@vuepress/plugin-back-to-top",
    "@vuepress/plugin-medium-zoom",
    ["vuepress-plugin-mathjax", { target: "svg", macros: { "*": "\\times" } }],
    [
      "vuepress-plugin-code-copy",
      {
        color: "#3eaf7c",
        backgroundTransition: false,
        staticIcon: false,
      },
    ],
    [
      '@vuepress/google-analytics',
      {
        'ga': 'UA-181048107-1' // the Google Analytics Tracking ID
      }
    ]
  ],
};
