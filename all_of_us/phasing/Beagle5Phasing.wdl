# WDL for phasing reference for Aou Local Ancestry
version 1.0

workflow Beagle5Phasing {
  input {
    String vcf_gs
    String map_gs
    String out_prefix
    String out_dir_gs

    Int window_markers = 3500000
    String java_xmx = "80g"

    File run_beagle_sh
  }
  String pipeline_version = "aou_9.0.0"

  call PhaseWithBeagle {
    input:
      vcf_gs = vcf_gs,
      map_gs = map_gs,
      out_prefix = out_prefix,
      out_dir_gs = out_dir_gs,
      window_markers = window_markers,
      java_xmx = java_xmx,
      run_beagle_sh = run_beagle_sh
  }

  output {
    String phased_vcf_gs = PhaseWithBeagle.phased_vcf_gs
    String phased_log_gs = PhaseWithBeagle.phased_log_gs
  }
}

task PhaseWithBeagle {
  input {
    String vcf_gs
    String map_gs
    String out_prefix
    String out_dir_gs

    Int window_markers = 3500000
    String java_xmx
	  File run_beagle_sh
	  Int runtime_cpu = 32
    Int mem_gb = 120
	  Int runtime_disk_gb = 2000
    String runtime_disk_type = "HDD"   # or "SSD"
  }

  command <<<
    chmod +x "~{run_beagle_sh}"
    bash beagle5_phasing.sh \
        --vcf "~{vcf_gs}" \
        --map "~{map_gs}" \
        --out-prefix "~{out_prefix}" \
        --out-dir "~{out_dir_gs}" \
        --window-markers "~{window_markers}" \
        --java-xmx "~{java_xmx}"
  >>>
  output {
	String phased_vcf_gs = "~{out_dir_gs}/~{out_prefix}.phased.vcf.gz"
	String phased_log_gs = "~{out_dir_gs}/~{out_prefix}.phased.log"
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/beagle5:0.0.1"
	cpu: runtime_cpu
	memory: mem_gb + " GB"
	disks: "local-disk ~{runtime_disk_gb} ~{runtime_disk_type}"
    # plus whatever backend-specific runtime fields you use
  }
}
