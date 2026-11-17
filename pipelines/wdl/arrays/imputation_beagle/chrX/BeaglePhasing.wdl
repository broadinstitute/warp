version 1.0

# One unit of phasing work: Beagle runs once per Region.
# Today there are up to three (PAR1, non-PAR, PAR2). To scatter later, the same
# struct can describe chunks: more entries with overlapping intervals.
struct Region {
    String name
    String interval    # chrom:start-end, 1-based inclusive
    File map
}

workflow BeaglePhasing {
    input {
        # Position-sorted, bgzipped, all-diploid VCF and its tabix index
        File input_vcf
        File input_vcf_index

        File beagle_jar = "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/jar/beagle.jar"
        File bref3_jar = "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/jar/bref3.jar"

        # Genetic maps for the separate regions
        File plink_map_par1 = "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/maps/plink.chrchrX_par1.GRCh38.map"
        File plink_map_nonpar = "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/maps/plink.chrchrX.GRCh38.map"
        File plink_map_par2 = "gs://fc-secure-6e8dde96-a72b-47ce-af63-92cadc8fb754/maps/plink.chrchrX_par2.GRCh38.map"

        # Region intervals (GRCh38 PAR boundaries)
        String par1_interval = "chrX:1-2781479"
        String nonpar_interval = "chrX:2781480-155701382"
        String par2_interval = "chrX:155701383-156040895"

        # Which regions to extract and phase. Set to false for regions the input
        # does not cover (e.g. a 5 Mb snippet); a selected region with no records
        # fails in ExtractRegions, before any Beagle VM starts.
        Boolean run_par1 = true
        Boolean run_nonpar = true
        Boolean run_par2 = true

        Int window_markers = 3200000

        # Copy compressed blocks instead of re-parsing every record. Requires
        # BGZF inputs with identical sample columns. Set to false if concat rejects them.
        Boolean naive_concat = true

        String chrom = "chrX"
        String prefix = "aou-anvil"

        Int preemptible_tries = 3
        # Beagle cannot resume, and large VMs are preempted often
        Int beagle_preemptible_tries = 0
    }

    # ---- Region list -------------------------------------------------------
    # Everything downstream works on `selected_regions`, in genomic order.
    # SCATTER LATER: build this list from chunks instead (each non-PAR chunk's
    # interval including its overlap, all using plink_map_nonpar).
    Array[Region] regions = [
        {"name": "par1", "interval": par1_interval, "map": plink_map_par1},
        {"name": "nonpar", "interval": nonpar_interval, "map": plink_map_nonpar},
        {"name": "par2", "interval": par2_interval, "map": plink_map_par2}
    ]
    Array[Boolean] run_region = [run_par1, run_nonpar, run_par2]

    scatter (i in range(length(regions))) {
        if (run_region[i]) {
            Region selected_or_none = regions[i]
        }
    }
    Array[Region] selected_regions = select_all(selected_or_none)

    scatter (r in selected_regions) {
        String selected_names = r.name
        String selected_intervals = r.interval
    }

    # ---- Extract every selected region from one copy of the input ----------
    call ExtractRegions {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            region_names = selected_names,
            region_intervals = selected_intervals,
            prefix = prefix,
            preemptibles = preemptible_tries
    }

    # ---- Phase each region independently -----------------------------------
    scatter (i in range(length(selected_regions))) {
        call RunBeagle {
            input:
                input_vcf = ExtractRegions.region_vcfs[i],
                beagle_jar = beagle_jar,
                plink_map = selected_regions[i].map,
                prefix = prefix + "." + selected_regions[i].name,
                window_markers = window_markers,
                preemptibles = beagle_preemptible_tries
        }
    }

    # ---- Join ---------------------------------------------------------------
    # The regions do not overlap, so plain concatenation is correct.
    # SCATTER LATER: overlapping chunks need a ligation step here instead
    # (match each sample's phase across the overlap, trim, then concatenate).
    call ConcatVcfs {
        input:
            phased_vcfs = RunBeagle.phased_vcf,
            out_prefix = prefix + "." + chrom + ".merged",
            naive = naive_concat,
            preemptibles = preemptible_tries
    }

    # The three tasks below run in parallel on the merged VCF

    call RunBref3 {
        input:
            phased_vcf = ConcatVcfs.merged_vcf,
            bref3_jar = bref3_jar,
            chrom = chrom,
            prefix = prefix,
            preemptibles = preemptible_tries
    }

    call ExtractMarkers {
        input:
            phased_vcf = ConcatVcfs.merged_vcf,
            chrom = chrom,
            prefix = prefix,
            preemptibles = preemptible_tries
    }

    call IndexVcf {
        input:
            phased_vcf = ConcatVcfs.merged_vcf,
            preemptibles = preemptible_tries
    }

    output {
        File phased_vcf = ConcatVcfs.merged_vcf
        File phased_vcf_tbi = IndexVcf.vcf_tbi
        Array[File] phased_logs = RunBeagle.phased_log
        Array[File] beagle_gc_logs = RunBeagle.gc_log
        File bref3_out = RunBref3.bref3_file
        File markers_out = ExtractMarkers.markers_file
    }
}

task ExtractRegions {
    input {
        File input_vcf
        File input_vcf_index
        Array[String] region_names
        Array[String] region_intervals
        String prefix

        String machine_type = "c3d-standard-16"
        Int memory_gb = 64
        Int cpu = 16
        # Input plus the extracted copies (together at most about the input size)
        Int disk_size_gb = ceil(size(input_vcf, "GB") * 2.5) + 50
        String disk_type = "SSD"
        Int preemptibles
    }

    command <<<
        set -euo pipefail

        NAMES=(~{sep=" " region_names})
        INTERVALS=(~{sep=" " region_intervals})
        n=${#NAMES[@]}
        if [ "$n" -eq 0 ]; then
            echo "No regions selected: set at least one run flag to true" >&2
            exit 1
        fi

        # tabix needs the index next to the VCF, named <vcf>.tbi
        ln -s ~{input_vcf} input.vcf.bgz
        ln -s ~{input_vcf_index} input.vcf.bgz.tbi

        # Share the cores between the regions, which are extracted in parallel
        threads=$(( ~{cpu} / n ))
        if [ "$threads" -lt 1 ]; then threads=1; fi
        tabix_threads=""
        usage=$(tabix --help 2>&1 || true)
        case "$usage" in
            *--threads*) tabix_threads="--threads $threads" ;;
        esac

        extract_region() {
            local out=$1 interval=$2
            local range=${interval#*:}
            # tabix jumps to the region through the index and copies whole lines;
            # nothing parses the genotype columns. tabix returns records that
            # OVERLAP the interval, so an indel starting just before it would
            # land in two regions: keep only records whose POS is inside.
            # shellcheck disable=SC2086
            tabix -h $tabix_threads input.vcf.bgz "$interval" \
                | START=${range%-*} END=${range#*-} perl -ne '
                    if (/^#/) { print; next; }
                    my ($pos) = /^[^\t]*\t(\d+)/;
                    print if $pos >= $ENV{START} && $pos <= $ENV{END};
                ' \
                | bgzip -@ "$threads" > "$out"
        }

        pids=()
        outs=()
        for i in "${!NAMES[@]}"; do
            # Zero-padded index keeps the outputs in region order
            out=$(printf 'region_%03d.%s.%s.vcf.gz' "$i" "~{prefix}" "${NAMES[$i]}")
            echo "Extracting ${NAMES[$i]} (${INTERVALS[$i]}) -> $out"
            extract_region "$out" "${INTERVALS[$i]}" &
            pids+=($!)
            outs+=("$out")
        done
        for pid in "${pids[@]}"; do
            wait "$pid"
        done

        # Fail here, cheaply, if a selected region has no records
        empty=0
        for i in "${!NAMES[@]}"; do
            found=$( { bgzip -dc "${outs[$i]}" 2>/dev/null || true; } \
                | { grep -v -m 1 '^#' || true; } | head -c 1 | wc -c )
            if [ "$found" -eq 0 ]; then
                echo "Region ${NAMES[$i]} (${INTERVALS[$i]}) has no records; set its run flag to false" >&2
                empty=1
            fi
        done
        if [ "$empty" -ne 0 ]; then exit 1; fi
        ls -l region_*.vcf.gz
    >>>

    output {
        # glob returns files in name order, which the zero-padded index makes region order
        Array[File] region_vcfs = glob("region_*.vcf.gz")
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        # Overrides cpu and memory below; keep those inputs in sync with the machine type
        predefinedMachineType: machine_type
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " " + disk_type
        preemptible: preemptibles
    }
}

task RunBeagle {
    input {
        File input_vcf
        File beagle_jar
        File plink_map
        String prefix
        Int window_markers

        String machine_type = "c2d-highmem-112"
        Int memory_gb = 896
        Int cpu = 112
        Int disk_size_gb = ceil(size(input_vcf, "GB") * 3) + 50
        String disk_type = "SSD"
        # Heap as a percentage of the VM's actual RAM. The remainder covers JVM
        # off-heap memory (GC structures, thread stacks, buffers) and the OS.
        Int java_heap_percent = 85
        Int preemptibles
    }

    command <<<
        set -euo pipefail

        echo "Running Beagle on ~{prefix}..."
        free -g || true
        # The GC log records heap use after each collection, which shows how
        # much memory this region really needs (for right-sizing the machine).
        java -ea -XX:MaxRAMPercentage=~{java_heap_percent} \
            -Xlog:gc:file=~{prefix}.gc.log \
            -jar ~{beagle_jar} \
            gt=~{input_vcf} \
            map=~{plink_map} \
            out=~{prefix} \
            nthreads=~{cpu} \
            window-markers=~{window_markers}
    >>>

    output {
        File phased_vcf = "~{prefix}.vcf.gz"
        File phased_log = "~{prefix}.log"
        File gc_log = "~{prefix}.gc.log"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-gcr-public/base/jre:17-debian"
        # Overrides cpu and memory below; keep those inputs in sync with the machine type
        predefinedMachineType: machine_type
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " " + disk_type
        preemptible: preemptibles
    }
}

task ConcatVcfs {
    input {
        Array[File] phased_vcfs
        String out_prefix
        Boolean naive = true

        # Smallest C3D type with local SSD (8 vCPU, 32 GB). Its local SSD capacity is fixed.
        String machine_type = "c3d-standard-8-lssd"
        Int memory_gb = 32
        Int cpu = 8
        Int disk_size_gb = ceil(size(phased_vcfs, "GB") * 2) + 50
        String disk_type = "LOCAL"
        Int preemptibles
    }

    command <<<
        set -euo pipefail

        if [ "~{naive}" = "true" ]; then
            echo "Concatenating phased regions (naive block copy)..."
            bcftools concat --naive -o ~{out_prefix}.vcf.gz ~{sep=" " phased_vcfs}
        else
            echo "Concatenating phased regions (full parse)..."
            bcftools concat --threads ~{cpu} -O z -o ~{out_prefix}.vcf.gz ~{sep=" " phased_vcfs}
        fi
    >>>

    output {
        File merged_vcf = "~{out_prefix}.vcf.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        # Overrides cpu and memory below; keep those inputs in sync with the machine type
        predefinedMachineType: machine_type
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " " + disk_type
        preemptible: preemptibles
    }
}

task RunBref3 {
    input {
        File phased_vcf
        File bref3_jar
        String chrom
        String prefix

        String machine_type = "c3d-highmem-8"
        Int memory_gb = 64
        Int cpu = 8
        Int disk_size_gb = ceil(size(phased_vcf, "GB") * 2) + 20
        String disk_type = "SSD"
        Int java_heap_percent = 85
        Int preemptibles
    }

    command <<<
        set -euo pipefail

        echo "Generating bref3 for ~{chrom}..."
        java -ea -XX:MaxRAMPercentage=~{java_heap_percent} -jar ~{bref3_jar} ~{phased_vcf} > ~{chrom}.~{prefix}.bref3
    >>>

    output {
        File bref3_file = "~{chrom}.~{prefix}.bref3"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-gcr-public/base/jre:17-debian"
        # Overrides cpu and memory below; keep those inputs in sync with the machine type
        predefinedMachineType: machine_type
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " " + disk_type
        preemptible: preemptibles
    }
}

task ExtractMarkers {
    input {
        File phased_vcf
        String chrom
        String prefix

        # Smallest C3D type with local SSD (8 vCPU, 32 GB). Its local SSD capacity is fixed.
        String machine_type = "c3d-standard-8-lssd"
        Int memory_gb = 32
        Int cpu = 8
        Int disk_size_gb = ceil(size(phased_vcf, "GB") * 2) + 20
        String disk_type = "LOCAL"
        Int preemptibles
    }

    command <<<
        set -euo pipefail

        echo "Extracting markers for ~{chrom}..."
        # Multi-threaded decompression, and cut before grep so grep only sees
        # the first 9 columns rather than every multi-megabyte line.
        bgzip -dc -@ ~{cpu} ~{phased_vcf} \
            | cut -f1-9 \
            | grep -v '^##contig' \
            | bgzip -@ ~{cpu} > ~{chrom}.~{prefix}.markers.gz
    >>>

    output {
        File markers_file = "~{chrom}.~{prefix}.markers.gz"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        # Overrides cpu and memory below; keep those inputs in sync with the machine type
        predefinedMachineType: machine_type
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " " + disk_type
        preemptible: preemptibles
    }
}

task IndexVcf {
    input {
        File phased_vcf

        # Smallest C3D type with local SSD (8 vCPU, 32 GB). Its local SSD capacity is fixed.
        String machine_type = "c3d-standard-8-lssd"
        Int memory_gb = 32
        Int cpu = 8
        Int disk_size_gb = ceil(size(phased_vcf, "GB") * 2) + 20
        String disk_type = "LOCAL"
        Int preemptibles
    }

    command <<<
        set -euo pipefail

        ln -s ~{phased_vcf} input.vcf.gz

        echo "Indexing VCF using ~{cpu} threads..."
        bcftools index --tbi --threads ~{cpu - 1} input.vcf.gz
    >>>

    output {
        File vcf_tbi = "input.vcf.gz.tbi"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        # Overrides cpu and memory below; keep those inputs in sync with the machine type
        predefinedMachineType: machine_type
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_size_gb + " " + disk_type
        preemptible: preemptibles
    }
}
