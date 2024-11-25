import subprocess

with open('/nobackup/h_cqs/shengq2/program/warp/.dockstore.yml', "rt") as fin:
  for line in fin:
    if "primaryDescriptorPath:" in line:
      filepath = '/nobackup/h_cqs/shengq2/program/warp' + line.split(":")[1].strip()
      print(filepath)
      subprocess.run(['java', '-jar', '/data/cqs/softwares/wdl/womtool-84.jar', 'validate', filepath])
