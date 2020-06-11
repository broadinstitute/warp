version 1.0

task checkOptimusInput {
  input {
    String chemistry
    String counting_mode
    Boolean force_no_check
  }  

  meta {
    description: "checks optimus input values and fails the pipeline immediately"
  }

  command {
    set -e
    
    ## Set pass to true
    pass="true"

    ## Perform checks
    if [[ ! ("${chemistry}" == "tenX_v2" || "${chemistry}" == "tenX_v3") ]]
    then
      pass="false"
      echo "ERROR: Invalid value \"${chemistry}\" for input \"chemistry\""
    fi

    if [[ ! ("${counting_mode}" == "sc_rna" || "${counting_mode}" == "sn_rna") ]]
    then
      pass="false"
      echo "ERROR: Invalid value \"${counting_mode}\" for input \"counting_mode\""
    fi

    if [[ ${force_no_check} == "true" ]]
    then
       echo "force_no_check is set: Ignoring input checks"
       exit 0;
    fi

    ## fail if any tests failed, ignore if force_no_check is set
    if [[ $pass == "true" ]]
    then
      exit 0;
    else
      exit 1;
    fi
  }

  runtime {
    docker: "ubuntu:18.04"
    cpu: 1
    memory: "1 GiB"
    disks: "local-disk 1 HDD"
  }
  
}
