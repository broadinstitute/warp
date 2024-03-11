version 1.0

workflow VUMCBedInfo {
  input {
    File input_fam
    File input_bim
  }

  call BedInfo {
    input:
      input_fam = input_fam,
      input_bim = input_bim
  }

  output {
    Int num_samples = BedInfo.num_samples
    Int num_variants = BedInfo.num_variants
  }
}

task BedInfo {
  input {
    File input_fam
    File input_bim
  }

  Int disk_size = ceil(size(input_fam, "GB") + size(input_bim, "GB")) + 2

  command <<<

wc -l ~{input_fam} | cut -d ' ' -f 1 > num_samples.txt

wc -l ~{input_bim} | cut -d ' ' -f 1 > num_variants.txt

>>>

  runtime {
    docker: "ubuntu:20.04"
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GiB"
  }
  output {
    Int num_samples = read_int("num_samples.txt")
    Int num_variants = read_int("num_variants.txt")
  }
}
