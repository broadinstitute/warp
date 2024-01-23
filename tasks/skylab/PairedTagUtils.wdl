version 1.0

task PairedTagDemultiplex {
    input {
        File read1_fastq
        File read3_fastq
        File barcodes_fastq
        String input_id

        # using the latest build of upstools docker in GCR
        String docker = "us.gcr.io/broad-gotc-prod/upstools:1.0.0-2023.03.03-1703173526"

        # Runtime attributes
        Int mem_size = 8
        Int cpu = 1
        # TODO decided cpu
        # estimate that bam is approximately equal in size to fastq, add 20% buffer
        Int disk_size = ceil(2 * ( size(read1_fastq, "GiB") + size(read3_fastq, "GiB") + size(barcodes_fastq, "GiB") )) + 400
        Int preemptible = 3
    }

    meta {
        description: "Demultiplexes paired-tag ATAC fastq files that have a 3 bp preindex and adds the index to readnames."
    }

    parameter_meta {
        read1_fastq: "read 1 FASTQ files of paired reads -- forward reads"
        read3_fastq: "read 3 FASTQ files of paired reads -- reverse reads"
        barcodes_fastq: "read 2 FASTQ files which contains the cellular barcodes"
        docker: "(optional) the docker image containing the runtime environment for this task"
        mem_size: "(optional) the amount of memory (MiB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
        disk_size: "(optional) the amount of disk space (GiB) to provision for this task"
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    }

    command <<<

        set -e
        echo ~{read1_fastq}
        echo ~{barcodes_fastq}
        echo ~{read3_fastq}
        echo Renaming files
        mv ~{read1_fastq} "~{input_id}_R1.fq.gz"
        mv ~{barcodes_fastq} "~{input_id}_R2.fq.gz"
        mv ~{read3_fastq} "~{input_id}_R3.fq.gz"

        echo Running UPStools
        upstools sepType_DPT ~{input_id} 3
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_size} GiB"
        disks: "local-disk ${disk_size} HDD"
        preemptible: preemptible
    }

    output {
        File fastq1 = "~{input_id}_R1_prefix.fq.gz"
        File barcodes = "~{input_id}_R2_prefix.fq.gz"
        File fastq3 = "~{input_id}_R3_prefix.fq.gz"
    }
}

task AddBBTag {
    input {
        File bam
        String input_id

        # using the latest build of upstools docker in GCR
        String docker = "us.gcr.io/broad-gotc-prod/upstools:1.0.0-2023.03.03-1704300311"

        # Runtime attributes
        Int mem_size = 8
        Int cpu = 1
        # TODO decided cpu
        # estimate that bam is approximately equal in size to fastq, add 20% buffer
        Int disk_size = ceil(2 * ( size(bam, "GiB"))) + 100
        Int preemptible = 3
    }

    meta {
        description: "Demultiplexes paired-tag ATAC fastq files that have a 3 bp preindex and adds the index to readnames."
    }

    parameter_meta {
        bam: "BAM with aligned reads and barcode in the CB tag"
        input_id: "input ID"
        docker: "(optional) the docker image containing the runtime environment for this task"
        mem_size: "(optional) the amount of memory (MiB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
        disk_size: "(optional) the amount of disk space (GiB) to provision for this task"
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    }

    command <<<

        set -e
        echo "BAM file name is:"
        echo ~{bam}
        echo moving BAM
        mv ~{bam} ./~{input_id}.bam
        echo Running UPStools
        python3 /upstools/pyscripts/scifi.preindex_CB_to_BB.py --in ~{input_id}.bam
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_size} GiB"
        disks: "local-disk ${disk_size} HDD"
        preemptible: preemptible
    }

    output {
        File bb_bam = "~{input_id}.bam.BB.bam"
    }
}
task ParseBarcodes {
    input {
        File atac_h5ad
        File atac_fragment
        Int nthreads = 1
        String cpuPlatform = "Intel Cascade Lake"
    }

    String atac_base_name = basename(atac_h5ad, ".h5ad")
    String atac_fragment_base = basename(atac_fragment, ".tsv")

    Int machine_mem_mb = ceil((size(atac_h5ad, "MiB") + size(gex_h5ad, "MiB") + size(atac_fragment, "MiB")) * 3) + 10000
    Int disk =  ceil((size(atac_h5ad, "GiB") + size(gex_h5ad, "GiB") + size(atac_fragment, "GiB")) * 5) + 10

  parameter_meta {
      atac_h5ad: "The resulting h5ad from the ATAC workflow."
      atac_fragment: "The resulting fragment TSV from the ATAC workflow."
  }

  command <<<
      set -e pipefail

      python3 <<CODE

      # set parameters
      atac_h5ad = "~{atac_h5ad}"
      atac_fragment = "~{atac_fragment}"

    # import anndata to manipulate h5ad files
      import anndata as ad
      import pandas as pd
      print("Reading ATAC h5ad:")
      print("~{atac_h5ad}")
      print("Read ATAC fragment file:")
      print("~{atac_fragment}")
    
      atac_data = ad.read_h5ad("~{atac_h5ad}")
      atac_tsv = pd.read_csv("~{atac_fragment}", sep="\t", names=['chr','start', 'stop', 'barcode','n_reads'])
      whitelist_atac = pd.read_csv("~{atac_whitelist}", header=None, names=["atac_barcodes"])

      # Separate out CB and preindex in the h5ad
      print("Creating CB and preindex columns")
      for x in range(len(atac_data.obs.index)):
        CB = atac_data.obs.index[x][3:]
        preindex = atac_data.obs.index[x][:3]
        atac_data.obs.loc[atac_data.obs.index[x], 'CB'] = CB
        atac_data.obs.loc[atac_data.obs.index[x], 'preindex'] = preindex
      
      # Indentify sample barcodes assigned to more than one cell barcode
      print("Identifying sample barcodes asigned to multiple cell barcodes in h5ad")
      list = []
      for preindex in atac_data.obs.preindex:
        if preindex not in list:
          list.append(preindex)
          CB = atac_data.obs.loc[atac_data.obs['preindex'] == preindex, "CB"]
          if len(CB)>1:
            atac_data.obs.loc[atac_data.obs['preindex'] == preindex, "Duplicates"] = '1'
          else:
            atac_data.obs.loc[atac_data.obs['preindex'] == preindex, "Duplicates"] = '0'

      # Idenitfy the barcodes in the whitelist that match barcodes in datasets
      atac_data.write_h5ad("~{atac_base_name}.h5ad")
      df_fragment.to_csv("~{atac_fragment_base}.tsv", sep='\t', index=False, header = False)
      CODE
      
      # sorting the file
      echo "Sorting file"
      sort -k1,1V -k2,2n "~{atac_fragment_base}.tsv" > "~{atac_fragment_base}.sorted.tsv"
      echo "Starting bgzip"
      bgzip "~{atac_fragment_base}.sorted.tsv"
      echo "Starting tabix"
      tabix -s 1 -b 2 -e 3 "~{atac_fragment_base}.sorted.tsv.gz"

  >>>

  runtime {
      docker: "us.gcr.io/broad-gotc-prod/snapatac2:1.0.4-2.3.1-1700590229"
      disks: "local-disk ~{disk} HDD"
      memory: "${machine_mem_mb} MiB"
      cpu: nthreads
  }

  output {
      File atac_h5ad_file = "~{atac_base_name}.h5ad"
      File atac_fragment_tsv = "~{atac_fragment_base}.sorted.tsv.gz"
      File atac_fragment_tsv_tbi = "~{atac_fragment_base}.sorted.tsv.gz.tbi"
  }
}
