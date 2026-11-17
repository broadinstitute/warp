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
        # Position-sorted, bgzipped, all-diploid, single-chromosome VCF.
        # No index is needed: regions are located by binary search.
        File input_vcf

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
    # Compressed blocks are copied as-is; nothing is decompressed except the
    # header and the two blocks at each region's edges.
    call ExtractRegions {
        input:
            input_vcf = input_vcf,
            region_names = selected_names,
            region_intervals = selected_intervals,
            prefix = prefix,
            preemptibles = preemptible_tries
    }

    # ---- Index the extracted region VCFs -----------------------------------
    # Nothing downstream reads these indices, so this runs alongside Beagle
    # and never delays it.
    scatter (region_vcf in ExtractRegions.region_vcfs) {
        call IndexVcf as IndexRegionVcf {
            input:
                vcf = region_vcf,
                preemptibles = preemptible_tries
        }
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
            vcf = ConcatVcfs.merged_vcf,
            preemptibles = preemptible_tries
    }

    output {
        File phased_vcf = ConcatVcfs.merged_vcf
        File phased_vcf_tbi = IndexVcf.vcf_tbi
        # Unphased per-region inputs to Beagle and their indices, in region order
        Array[File] region_vcfs = ExtractRegions.region_vcfs
        Array[File] region_vcf_tbis = IndexRegionVcf.vcf_tbi
        Array[File] phased_logs = RunBeagle.phased_log
        Array[File] beagle_gc_logs = RunBeagle.gc_log
        File bref3_out = RunBref3.bref3_file
        File markers_out = ExtractMarkers.markers_file
    }
}

task ExtractRegions {
    input {
        File input_vcf
        Array[String] region_names
        Array[String] region_intervals
        String prefix

        # Almost all of the work is copying bytes, so a small machine is enough
        String machine_type = "c3d-standard-4"
        Int memory_gb = 16
        Int cpu = 4
        # Input plus the extracted copies (together at most about the input size)
        Int disk_size_gb = ceil(size(input_vcf, "GB") * 2.5) + 50
        String disk_type = "SSD"
        Int preemptibles
    }

    command <<<
        set -euo pipefail

        NAMES=(~{sep=" " region_names})
        INTERVALS=(~{sep=" " region_intervals})
        if [ "${#NAMES[@]}" -eq 0 ]; then
            echo "No regions selected: set at least one run flag to true" >&2
            exit 1
        fi

        cat > bgzf_extract.py << 'PYEOF'
        #!/usr/bin/env python3
        # Extract position intervals from a position-sorted, single-chromosome BGZF VCF
        # by copying compressed blocks. No index is used. For each interval:
        #   1. binary-search the compressed file for the first record at or after the
        #      interval start, and the first record after the interval end
        #   2. write a fresh header, then copy every whole block between the two points
        #      byte for byte; only the two boundary blocks are decompressed and re-cut
        # Usage: bgzf_extract.py IN.vcf.bgz OUT1 chrom:start-end [OUT2 chrom:start-end ...]
        import struct
        import sys
        import zlib

        NL = bytes([10])
        TAB = bytes([9])
        MAGIC = bytes.fromhex("1f8b08040000000000ff060042430200")
        EOF_BLOCK = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
        HEAD = 64            # bytes of each line needed to read CHROM and POS
        WINDOW = 256 * 1024  # stop bisecting when the search window is this small
        COPY = 16 * 1024 * 1024


        def read_block(f, off):
            # Return (data, next_off) for the block at off, or None if there is none.
            f.seek(off)
            hdr = f.read(18)
            if len(hdr) < 18 or hdr[:16] != MAGIC:
                return None
            size = struct.unpack("<H", hdr[16:18])[0] + 1
            raw = hdr + f.read(size - 18)
            if len(raw) < size:
                return None
            try:
                return zlib.decompress(raw, 31), off + size
            except zlib.error:
                return None


        def find_block(f, off, file_size):
            # Offset of the first real block starting at or after off (file_size if none).
            while off < file_size:
                f.seek(off)
                buf = f.read(WINDOW + 16)
                i = buf.find(MAGIC)
                while i != -1:
                    if read_block(f, off + i) is not None:
                        return off + i
                    i = buf.find(MAGIC, i + 1)
                if len(buf) < WINDOW + 16:
                    break
                off += WINDOW
            return file_size


        def line_starts(f, off, file_size):
            # Yield (block_off, intra, head) for every line that starts at or after the
            # first newline found from block offset off. head holds the line's first bytes.
            pending = None
            off = find_block(f, off, file_size)
            while off < file_size:
                blk = read_block(f, off)
                if blk is None:
                    raise SystemExit("corrupt BGZF block at offset %d" % off)
                data, next_off = blk
                if pending is not None:
                    pending[2] += data[:HEAD - len(pending[2])]
                    if len(pending[2]) >= HEAD or NL in pending[2]:
                        yield pending[0], pending[1], bytes(pending[2])
                        pending = None
                i = 0
                while True:
                    j = data.find(NL, i)
                    if j == -1:
                        break
                    s = j + 1
                    if s == len(data):
                        pending = [next_off, 0, bytearray()]
                    else:
                        head = data[s:s + HEAD]
                        if len(head) >= HEAD or NL in head:
                            yield off, s, head
                        else:
                            pending = [off, s, bytearray(head)]
                    i = s
                off = next_off


        def line_pos(head, chrom):
            # -1 for header lines, else POS. The file must hold a single chromosome.
            if head[:1] == b"#":
                return -1
            fields = head.split(TAB)
            if len(fields) < 3:
                raise SystemExit("cannot read POS from line starting %r" % head[:40])
            if fields[0].decode() != chrom:
                raise SystemExit("found contig %s, expected only %s" % (fields[0].decode(), chrom))
            return int(fields[1])


        def split_point(f, file_size, chrom, boundary):
            # (block_off, intra, pos) of the first record with POS >= boundary,
            # or (file_size, 0, None) if there is none.
            lo, hi = 0, file_size
            while hi - lo > WINDOW:
                mid = (lo + hi) // 2
                first = next(line_starts(f, mid, file_size), None)
                if first is None or line_pos(first[2], chrom) >= boundary:
                    hi = mid
                else:
                    lo = mid
            for block_off, intra, head in line_starts(f, lo, file_size):
                pos = line_pos(head, chrom)
                if pos >= boundary:
                    return block_off, intra, pos
            return file_size, 0, None


        def write_blocks(out, data):
            for i in range(0, len(data), 60000):
                chunk = data[i:i + 60000]
                comp = zlib.compressobj(6, zlib.DEFLATED, -15)
                body = comp.compress(chunk) + comp.flush()
                out.write(MAGIC + struct.pack("<H", len(body) + 25) + body
                          + struct.pack("<II", zlib.crc32(chunk), len(chunk)))


        def copy_raw(f, out, start, end):
            f.seek(start)
            left = end - start
            while left > 0:
                buf = f.read(min(COPY, left))
                if not buf:
                    raise SystemExit("unexpected end of file while copying")
                out.write(buf)
                left -= len(buf)


        def main():
            path = sys.argv[1]
            jobs = sys.argv[2:]
            if not jobs or len(jobs) % 2:
                raise SystemExit(__doc__ or "usage: IN OUT interval [OUT interval ...]")
            with open(path, "rb") as f:
                file_size = f.seek(0, 2)
                # A trailing empty EOF block is not data: never copy it
                data_end = file_size
                if file_size >= 28:
                    f.seek(file_size - 28)
                    if f.read(28) == EOF_BLOCK:
                        data_end = file_size - 28

                # Header = everything before the first record
                chrom0 = jobs[1].split(":")[0]
                h_block, h_intra, _ = split_point(f, data_end, chrom0, 0)
                header = bytearray()
                off = 0
                while off < h_block:
                    data, off = read_block(f, off)
                    header += data
                if h_intra:
                    header += read_block(f, h_block)[0][:h_intra]

                failed = False
                for k in range(0, len(jobs), 2):
                    out_path, interval = jobs[k], jobs[k + 1]
                    chrom, rng = interval.split(":")
                    start, end = (int(x) for x in rng.split("-"))
                    a_block, a_intra, a_pos = split_point(f, data_end, chrom, start)
                    b_block, b_intra, _ = split_point(f, data_end, chrom, end + 1)
                    if (a_block, a_intra) == (b_block, b_intra):
                        print("%s: no records in %s" % (out_path, interval), file=sys.stderr)
                        failed = True
                        continue
                    with open(out_path, "wb") as out:
                        write_blocks(out, bytes(header))
                        if a_block == b_block:
                            write_blocks(out, read_block(f, a_block)[0][a_intra:b_intra])
                        else:
                            raw_start = a_block
                            if a_intra:
                                data, raw_start = read_block(f, a_block)
                                write_blocks(out, data[a_intra:])
                            copy_raw(f, out, raw_start, b_block)
                            if b_intra:
                                write_blocks(out, read_block(f, b_block)[0][:b_intra])
                        out.write(EOF_BLOCK)
                    print("%s: %s, first record at %d, %d compressed bytes copied"
                          % (out_path, interval, a_pos, max(0, b_block - a_block)), file=sys.stderr)
                if failed:
                    raise SystemExit("a selected region has no records; set its run flag to false")


        if __name__ == "__main__":
            main()
        PYEOF

        args=()
        for i in "${!NAMES[@]}"; do
            # Zero-padded index keeps the outputs in region order
            out=$(printf 'region_%03d.%s.%s.vcf.gz' "$i" "~{prefix}" "${NAMES[$i]}")
            args+=("$out" "${INTERVALS[$i]}")
        done

        # Fails if a selected region has no records
        python3 bgzf_extract.py ~{input_vcf} "${args[@]}"
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

        String machine_type = "c3d-standard-4"
        Int memory_gb = 16
        Int cpu = 4
        Int disk_size_gb = ceil(size(phased_vcfs, "GB") * 2) + 50
        # Persistent SSD: the local SSD of the -lssd machine types was not
        # usable as the Cromwell work disk (tasks ran out of disk)
        String disk_type = "SSD"
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

        String machine_type = "c3d-standard-4"
        Int memory_gb = 16
        Int cpu = 4
        Int disk_size_gb = ceil(size(phased_vcf, "GB") * 2) + 20
        # Persistent SSD: the local SSD of the -lssd machine types was not
        # usable as the Cromwell work disk (tasks ran out of disk)
        String disk_type = "SSD"
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
        File vcf

        String machine_type = "c3d-standard-8"
        Int memory_gb = 32
        Int cpu = 8
        Int disk_size_gb = ceil(size(vcf, "GB") * 2) + 20
        # Persistent SSD: the local SSD of the -lssd machine types was not
        # usable as the Cromwell work disk (tasks ran out of disk)
        String disk_type = "SSD"
        Int preemptibles
    }

    # The index is named after the VCF, so it can sit next to it later
    String vcf_name = basename(vcf)

    command <<<
        set -euo pipefail

        ln -s ~{vcf} ~{vcf_name}

        echo "Indexing ~{vcf_name} using ~{cpu} threads..."
        bcftools index --tbi --threads ~{cpu - 1} ~{vcf_name}
    >>>

    output {
        File vcf_tbi = "~{vcf_name}.tbi"
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
