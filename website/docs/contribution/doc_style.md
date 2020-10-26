# Documentation Style Guide

This guide provides some examples about how to add new documentation that can be properly rendered on this website. Please note most of the Github flavored [Markdown](https://github.github.com/gfm/) syntax should work natrually, this guide just tries to elaboratethe extension syntax to it.

## 1. Insert code blocks

This site supports inserting code blocks with highlighted lines, for examples, the following:

````md
```wdl {3-4,7-11}
task SortBam {
    input {
        File bam_input
        String sort_order = "coordinate"

        # runtime values
        String docker = "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
        Int machine_mem_mb = 8250
        Int machine_overhead_mb = 500
        Int cpu = 1
        Int preemptible = 3
    }

    Int command_mem_mb = machine_mem_mb - machine_overhead_mb
    Int disk = ceil(size(bam_input, "Gi") * 6) + 50

    meta {
        description: "Sorts bam"
    }

    command {
        set -e

        java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar SortSam \
              I=${bam_input} \
              O=sorted.bam \
              SORT_ORDER=${sort_order}
    }

    runtime {
        docker: docker
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File bam_output = "sorted.bam"
    }
}
```
````

will result in a rendered code block like:

```wdl {3-4,7-11}
task SortBam {
    input {
        File bam_input
        String sort_order = "coordinate"

        # runtime values
        String docker = "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
        Int machine_mem_mb = 8250
        Int machine_overhead_mb = 500
        Int cpu = 1
        Int preemptible = 3
    }

    Int command_mem_mb = machine_mem_mb - machine_overhead_mb
    Int disk = ceil(size(bam_input, "Gi") * 6) + 50

    meta {
        description: "Sorts bam"
    }

    command {
        set -e

        java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar SortSam \
              I=${bam_input} \
              O=sorted.bam \
              SORT_ORDER=${sort_order}
    }

    runtime {
        docker: docker
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File bam_output = "sorted.bam"
    }
}
```

## 2. Insert tables

```md
| Some Table Col 1 | Some Table Col 2 |
| :--------------: | :--------------: |
|       Val1       |       Val4       |
|       Val2       |       Val5       |
|       Val3       |       Val6       |
```

| Some Table Col 1 | Some Table Col 2 |
| :--------------: | :--------------: |
|       Val1       |       Val4       |
|       Val2       |       Val5       |
|       Val3       |       Val6       |

::: tip TIP
It's worth mentioning that [Tables Generator](https://www.tablesgenerator.com/markdown_tables) is a great tool for generating and re-formatting markdown tables.
:::

## 3. Cross-reference and anchor

To link to another section within the same article, you would use `[Return to ## 1. Insert code blocks](#_1-insert-code-blocks)`: [Return to ## 1. Insert code blocks](#_1-insert-code-blocks).

To link to sections in other articles, use the following syntax (note the relative paths):

- `[Return to Changelog Style Guide](./changelog_style)`: [Return to Changelog Style Guide](./changelog_style)
- `[Return to The Documentation](../documentation/About_WARP/BestPractices.md#Best-Practices-for-Building-Data-Processing Pipelines)`: [Return to The Documentation](../documentation/About_WARP/BestPractices.md#Best-Practices-for-Building-Data-Processing-Pipelines)

## 4. Centered text block

To make a text block centered, use:

```md
<center>
Centered Text Block!
</center>
```

<center>
Centered Text Block!
</center>


::: danger NOTE
For images, you **HAVE TO** insert blank lines to make them work:

```md
<center>

![](./some_pic.png)

</center>
```
:::

## 5. Text with color backgorund

You could use the following to highlight your text:

```html
<span id="inline-blue"> Text with blue background </span>,
<span id="inline-purple"> Text with purple background </span>,
<span id="inline-yellow"> Text with yellow background </span>,
<span id="inline-green"> Text with green background </span>
```

<span id="inline-blue"> Text with blue background </span>,
<span id="inline-purple"> Text with purple background </span>,
<span id="inline-yellow"> Text with yellow background </span>,
<span id="inline-green"> Text with green background </span>

You could also use the markdown-it syntax `==My text==` to achieve similar result: ==My text==.

## 6. Custom containers

As we already saw in this guide several places, we could add custom containers:


```md
::: tip
This is a tip without title!
:::
```

::: tip
This is a tip without title!
:::

```md
::: tip
This is a tip with a title!
:::
```

::: tip TITLE
This is a tip with a title!
:::


```md WARNING
::: warning
This is a warning!
:::
```

::: warning WARNING
This is a warning!
:::

```md
::: danger DANGER
This is a danger!
:::
```

::: danger DANGER
This is a danger!
:::

```md
::: details SUMMARY
Details by default are hidden!
:::
```

::: details SUMMARY
Details by default are hidden!
:::

## 7. Code groups

You could also insert tab-based code groups:

<code-group>
<code-block title="Python">
```python
print("Hello World")
```
</code-block>

<code-block title="C++">
```cpp
#include<iostream>

using namespace std;

int main()
{
    // prints hello world
    cout<<"Hello World";

    return 0;
}
```
</code-block>

<code-block title="R">
```r
print("Hello World")
```
</code-block>
</code-group>

## 8. Footnotes

It is important to cite the references, to do so, use the `markdown-it`'s footnotes syntax:

```md
This sentence has a footnote[^1]. (See footnote at the bottom of this guide.)

[^1]: I'm a footnote!
```

which results in:

---

This sentence has a footnote[^1]. (See footnote at the bottom of this guide.)

[^1]: I'm a footnote!

---

We could also write in-line footnotes, which is much easier to write without counting back and forth:

```md
This sentence has another footnote ^[I'm another footnote] (See footnote at the bottom of this page.)
```

which has the same effect:

---

This sentence has another footnote ^[I'm another footnote] (See footnote at the bottom of this page.)

---

## 9. Insert images

Insert images is as straight-forward as using the ordinary markdown syntax:

```md
![terra](./Terra_warp.png)
```

![terra](./Terra_warp.png)

you could also use syntax like the following to resize:

```md
![terra](./Terra_warp.png =200x60)
```

![terra](./Terra_warp.png =200x60)

## 10. Insert Table of Contents (ToC)

You could use:

```md
[[toc]]
```

to insert in-line ToC:

[[toc]]

## 11. Insert badges

You could add a variety of badges along with the text:

```md
<Badge text="hello!" type="tip"/>
<Badge text="hello!" type="tip" vertical="middle"/>
```
<Badge text="hello!" type="tip"/>
<Badge text="hello!" type="tip" vertical="middle"/>

```md
<Badge text="hello!" type="warning"/>
<Badge text="hello!" type="warning" vertical="middle"/>
```
<Badge text="hello!" type="warning"/>
<Badge text="hello!" type="warning" vertical="middle"/>

```md
<Badge text="hello!" type="error"/>
<Badge text="hello!" type="error" vertical="middle"/>
```
<Badge text="hello!" type="error"/>
<Badge text="hello!" type="error" vertical="middle"/>
