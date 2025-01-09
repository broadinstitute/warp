version 1.0
workflow BuildIndexHisat3n {
	input {
		File genome_fa
		File genome_ver
		Int memory
        Int disk = ceil(size(genome_fa, "GiB"))*2 + 50
        String hisat_command
        File? monitoring_script_path
	}
	call hisat_build_methyl {
		input:
		  genome_fa = genome_fa,
		  genome_ver = genome_ver,
		  memory = memory,
          disk = disk,
          hisat_command = hisat_command
	}
	output {
		Array[File] final_out = hisat_build_methyl.output_files
        File monitoring_hisat3 = hisat_build_methyl.monitoring_out
	}
}

task hisat_build_methyl {
  input {
    File genome_fa
    String genome_ver
    Int memory
    Int disk
    String hisat_command = "/hisat-3n/hisat-3n-build --base-change C,T --repeat-index"
    File monitoring = "gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"
  }
  command <<<
    set -euo pipefail
    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring}" ]; then
      chmod a+x ~{monitoring}
      ~{monitoring} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi
    ~{hisat_command} ~{genome_fa} ~{genome_ver}
  >>>
  output {
    Array[File] output_files = glob("hg38*")
    File monitoring_out = "monitoring.log"
  }
  runtime {
    docker: "ekiernan/hisat3n-python:v1"
    memory: memory + " GiB"
    disks: "local-disk ~{disk} HDD"
    disk: disk + " GB" # TES
  }
}