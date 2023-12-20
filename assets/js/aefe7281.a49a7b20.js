"use strict";(self.webpackChunkwebsite_2=self.webpackChunkwebsite_2||[]).push([[6461],{3905:function(e,n,t){t.d(n,{Zo:function(){return h},kt:function(){return u}});var a=t(7294);function i(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function o(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);n&&(a=a.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,a)}return t}function r(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?o(Object(t),!0).forEach((function(n){i(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):o(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function l(e,n){if(null==e)return{};var t,a,i=function(e,n){if(null==e)return{};var t,a,i={},o=Object.keys(e);for(a=0;a<o.length;a++)t=o[a],n.indexOf(t)>=0||(i[t]=e[t]);return i}(e,n);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(a=0;a<o.length;a++)t=o[a],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(i[t]=e[t])}return i}var s=a.createContext({}),d=function(e){var n=a.useContext(s),t=n;return e&&(t="function"==typeof e?e(n):r(r({},n),e)),t},h=function(e){var n=d(e.components);return a.createElement(s.Provider,{value:n},e.children)},c="mdxType",p={inlineCode:"code",wrapper:function(e){var n=e.children;return a.createElement(a.Fragment,{},n)}},g=a.forwardRef((function(e,n){var t=e.components,i=e.mdxType,o=e.originalType,s=e.parentName,h=l(e,["components","mdxType","originalType","parentName"]),c=d(t),g=i,u=c["".concat(s,".").concat(g)]||c[g]||p[g]||o;return t?a.createElement(u,r(r({ref:n},h),{},{components:t})):a.createElement(u,r({ref:n},h))}));function u(e,n){var t=arguments,i=n&&n.mdxType;if("string"==typeof e||i){var o=t.length,r=new Array(o);r[0]=g;var l={};for(var s in n)hasOwnProperty.call(n,s)&&(l[s]=n[s]);l.originalType=e,l[c]="string"==typeof e?e:i,r[1]=l;for(var d=2;d<o;d++)r[d]=t[d];return a.createElement.apply(null,r)}return a.createElement.apply(null,t)}g.displayName="MDXCreateElement"},4832:function(e,n,t){t.r(n),t.d(n,{assets:function(){return s},contentTitle:function(){return r},default:function(){return p},frontMatter:function(){return o},metadata:function(){return l},toc:function(){return d}});var a=t(3117),i=(t(7294),t(3905));const o={sidebar_position:1},r="Changelog Style Guide Overview",l={unversionedId:"contribution/contribute_to_warp/changelog_style",id:"contribution/contribute_to_warp/changelog_style",title:"Changelog Style Guide Overview",description:"The style guide provides formatting guidelines and language suggestions for pipeline changelogs. It is divided into four sections: 1) Changelog Information, which details the types of changes listed in a changelog, 2) Language Usage, which describes language style and syntax for changelog information, 3) Changelog Format, which provides sample formatting for changelog information, and 4) Sample Changelog Entries, which displays two examples of changelog entries taken from the Optimus.changelog.md file.",source:"@site/docs/contribution/contribute_to_warp/changelog_style.md",sourceDirName:"contribution/contribute_to_warp",slug:"/contribution/contribute_to_warp/changelog_style",permalink:"/warp/docs/contribution/contribute_to_warp/changelog_style",draft:!1,editUrl:"https://github.com/broadinstitute/warp/edit/develop/website/docs/contribution/contribute_to_warp/changelog_style.md",tags:[],version:"current",lastUpdatedBy:"Kaylee Mathews",lastUpdatedAt:1703084290,formattedLastUpdatedAt:"Dec 20, 2023",sidebarPosition:1,frontMatter:{sidebar_position:1},sidebar:"docsSidebar",previous:{title:"Welcome",permalink:"/warp/docs/contribution/README"},next:{title:"WDL Task Runtime Style Guide Overview",permalink:"/warp/docs/contribution/contribute_to_warp/wdl_task_runtime_style"}},s={},d=[{value:"Changelog Information",id:"changelog-information",level:2},{value:"1) The changelog.md file name",id:"1-the-changelogmd-file-name",level:3},{value:"2) The pipeline release version name",id:"2-the-pipeline-release-version-name",level:3},{value:"3) The date of last commit (YYYY-MM-DD)",id:"3-the-date-of-last-commit-yyyy-mm-dd",level:3},{value:"4) Details of the version changes",id:"4-details-of-the-version-changes",level:3},{value:"Language Usage",id:"language-usage",level:2},{value:"Changelog Format",id:"changelog-format",level:2},{value:"Sample Changelog Entries",id:"sample-changelog-entries",level:2},{value:"Major Version Update (Breaking Change)",id:"major-version-update-breaking-change",level:3},{value:"Breaking changes",id:"breaking-changes",level:3},{value:"Additional changes",id:"additional-changes",level:3},{value:"Non-major Version Update (Non-breaking Change)",id:"non-major-version-update-non-breaking-change",level:3}],h={toc:d},c="wrapper";function p(e){let{components:n,...t}=e;return(0,i.kt)(c,(0,a.Z)({},h,t,{components:n,mdxType:"MDXLayout"}),(0,i.kt)("h1",{id:"changelog-style-guide-overview"},"Changelog Style Guide Overview"),(0,i.kt)("p",null,"The style guide provides formatting guidelines and language suggestions for pipeline changelogs. It is divided into four sections: 1) Changelog Information, which details the types of changes listed in a changelog, 2) Language Usage, which describes language style and syntax for changelog information, 3) Changelog Format, which provides sample formatting for changelog information, and 4) Sample Changelog Entries, which displays two examples of changelog entries taken from the ",(0,i.kt)("a",{parentName:"p",href:"https://github.com/broadinstitute/warp/blob/develop/pipelines/skylab/optimus/Optimus.changelog.md"},"Optimus.changelog.md file"),"."),(0,i.kt)("h2",{id:"changelog-information"},"Changelog Information"),(0,i.kt)("p",null,"Pipeline changelogs have four informational parts:"),(0,i.kt)("h3",{id:"1-the-changelogmd-file-name"},"1) The changelog.md file name"),(0,i.kt)("p",null,"The changelog file name should be consistent across pipelines. The format is ",(0,i.kt)("inlineCode",{parentName:"p"},"<pipeline name>.changelog.md"),"."),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"Ex: Optimus.changelog.md")),(0,i.kt)("p",null,"The file name is not included in the actual changelog- it is just the file name."),(0,i.kt)("h3",{id:"2-the-pipeline-release-version-name"},"2) The pipeline release version name"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"Ex: ",(0,i.kt)("a",{parentName:"li",href:"https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/genotyping/illumina/IlluminaGenotypingArray.changelog.md"},"Illumina Genotyping Array pipeline"),' version name is a number (i.e. "1.0.0")')),(0,i.kt)("p",null,"The version name is listed at the top of each changelog entry section. It should follow ",(0,i.kt)("a",{parentName:"p",href:"https://semver.org/"},"Semantic versioning (SemVer)")," which uses the major.minor.patch version number."),(0,i.kt)("h3",{id:"3-the-date-of-last-commit-yyyy-mm-dd"},"3) The date of last commit (YYYY-MM-DD)"),(0,i.kt)("p",null,"The purpose of the date of last commit is to inform users of the relative timing of changes to the pipeline. This is ",(0,i.kt)("strong",{parentName:"p"},"not a release date"),", which would have to be added to the changelog post-release. This date is updated as each change is created. For each pipeline release in the changelog, there will only be one final commit date (as opposed to a commit date for every change in that particular release)."),(0,i.kt)("h3",{id:"4-details-of-the-version-changes"},"4) Details of the version changes"),(0,i.kt)("p",null,"This section describes (in bullet point format) every type of change made to the pipeline for the current pipeline version. Types of changes include: additions, updates, changes, deprecations, removals, fixes, and security updates. The type of change should be used as the first word of the bullet point (see examples below). These changes should be updated as pipeline changes are made even if the updated pipeline version has not been released. This will enable end-users to see what changes are coming when the new release is published."),(0,i.kt)("p",null,'If changes are breaking changes to the pipeline (a major version update), this changelog section will be divided into two subsections: "Breaking changes" and "Additional changes". An example of these headers is shown in the ',(0,i.kt)("a",{parentName:"p",href:"#sample-changelog-entries"},"Sample Changelog Entries")," section."),(0,i.kt)("h2",{id:"language-usage"},"Language Usage"),(0,i.kt)("p",null,"All changelog notes should be bulleted (no period at the end of the final sentence of any one bullet point). Each bullet should include one type of change, but more than one sentence can be used to describe the individual change. Bullets should be written in complete sentences, but with the caveat that personal pronouns (\u201cI\u201d, \u201cWe\u201d, \u201cThey\u201d) are omitted. The first word of each bullet should be capitalized and describe the type of change made (in past tense)."),(0,i.kt)("p",null,"Examples of bullet points:"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},'"Updated the WDL with the latest docker image\u201d'),(0,i.kt)("li",{parentName:"ul"},'"Added an output to the WDL. Users will now see this output after running the WDL"')),(0,i.kt)("p",null,"For all bullet points, use active voice instead of passive voice. Examples are shown below."),(0,i.kt)("p",null,"Examples of passive voice:"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},'"Broken link in the StarAlign task ',(0,i.kt)("strong",{parentName:"li"},"was fixed"),'"')),(0,i.kt)("p",null,"Example of active voice:"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},'"',(0,i.kt)("strong",{parentName:"li"},"Fixed"),' the broken link in the StarAlign task"')),(0,i.kt)("p",null,"You can find more examples of active voice from the ",(0,i.kt)("a",{parentName:"p",href:"https://writing.wisc.edu/handbook/style/ccs_activevoice/"},"University of Wisconsin-Madison's Writing Center"),"."),(0,i.kt)("h2",{id:"changelog-format"},"Changelog Format"),(0,i.kt)("p",null,"The following is the markdown format required for all major changelog version updates which have breaking changes:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-md"},'# "insert Release Version Name here"\n\n"insert date of last commit in YYYY-MM-DD here" (Date of Last Commit)\n\n### Breaking changes\n\n* "insert description of first breaking pipeline change here"\n* "insert description of additional breaking changes here- add more bullets as necessary for additional changes"\n\n### Additional changes\n\n* "insert description of first non-breaking pipeline change here"\n* "insert description of additional change here- add more bullets as necessary for additional changes"\n')),(0,i.kt)("p",null,"The following is the markdown format required for non-major changelog version updates with non-breaking changes:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-md"},'# "insert Release Version Name here"\n\n"insert date of last commit in YYYY-MM-DD here" (Date of Last Commit)\n\n* "insert description of first pipeline change here"\n* "insert description of additional change here- add more bullets as necessary for additional changes"\n')),(0,i.kt)("h2",{id:"sample-changelog-entries"},"Sample Changelog Entries"),(0,i.kt)("h3",{id:"major-version-update-breaking-change"},"Major Version Update (Breaking Change)"),(0,i.kt)("h1",{id:"400"},"4.0.0"),(0,i.kt)("p",null,"2020-08-10 (Date of Last Commit)"),(0,i.kt)("h3",{id:"breaking-changes"},"Breaking changes"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"Changed sample_id to input_id")),(0,i.kt)("h3",{id:"additional-changes"},"Additional changes"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"Added input_name as an optional input for user provided sample_id"),(0,i.kt)("li",{parentName:"ul"},"Passed pipeline_version to output loom file"),(0,i.kt)("li",{parentName:"ul"},"Added input_id_metadata_field and input_name_metadata_field as optional input")),(0,i.kt)("h3",{id:"non-major-version-update-non-breaking-change"},"Non-major Version Update (Non-breaking Change)"),(0,i.kt)("h1",{id:"140"},"1.4.0"),(0,i.kt)("p",null,"2019-11-08 (Date of Last Commit)"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"Added support for V3 chemistry"),(0,i.kt)("li",{parentName:"ul"},"Updated the documentation with additional information for the README and additional files for Loom schema and BAM tags"),(0,i.kt)("li",{parentName:"ul"},"Updated the Zarr output")))}p.isMDXComponent=!0}}]);