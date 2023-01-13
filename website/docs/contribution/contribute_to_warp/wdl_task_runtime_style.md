---
sidebar_position: 2
---

# WDL Task Runtime Style Guide Overview

This style guide provides formatting guidelines and best practices for the runtime block of a WDL workflow task. For more information about scripting in WDL, see the [WDL 1.0 Specification](https://github.com/openwdl/wdl/blob/main/SPEC.md).

## Variables and suggested configurations
- Disk
    - Set in GiB
    - Variable type is Int
    - Use `local-disk` and `HDD`
    - Default should be dynamically set based on task inputs
    - Variable name should be `disk_size_gb`
- Memory
    - Set in MiB
    - Variable type is Int
    - Default should be dynamically set based on task inputs
    - Variable name should be `memory_mb`
- CPU
    - Default set to 1
    - Variable name should be `cpu`
- Docker
    - Variable name should be `docker`

## Example task input

Include runtime variables in the task inputs. Always provide a default value. Separate runtime variables from task inputs with one blank line.

```
    input {
        File example_file

        String docker = "us.gcr.io/path/to/docker"
        Int cpu = 1
        Int memory_mb = ceil(size(example_file, "MiB"))
        Int disk_size_db = ceil(size(example_file, "GiB"))
    }
```

## Example task runtime

Runtime block should be positioned between the command and output blocks of the task.

```
    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks: "local-disk ${disk_size_gb} HDD"
    }
```