version 1.0

task string_to_array {
  input {
    String str
    String delimiter = ","
  }
  command {
    echo ~{str} | tr '~{delimiter}' '\n'
  }
  output {
    Array[String] arr = read_lines(stdout())
  }
}
