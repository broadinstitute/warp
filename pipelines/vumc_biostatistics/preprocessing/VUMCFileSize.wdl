version 1.0

workflow VUMCFileSize {
  input {
    File input_file1
    File? input_file2
    File? input_file3
    File? input_file4
    File? input_file5
    File? input_file6
  }

  Float file1_size = size(input_file1)

  if(defined(input_file2)){
    Float file2_size =  size(input_file2)
  }

  if(defined(input_file3)){
    Float file3_size =  size(input_file3)
  }

  if(defined(input_file4)){
    Float file4_size =  size(input_file4)
  }

  if(defined(input_file5)){
    Float file5_size =  size(input_file5)
  }

  if(defined(input_file6)){
    Float file6_size =  size(input_file6)
  }

  output {
    Float input_file1_size = file1_size
    Float? input_file2_size = file2_size
    Float? input_file3_size = file3_size
    Float? input_file4_size = file4_size
    Float? input_file5_size = file5_size
    Float? input_file6_size = file6_size
  }
}
