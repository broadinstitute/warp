version 1.0 

workflow CaclulateAF {
    input {
        File pvar
        File psam
        File pgen
        File sample_list
        String prefix
        Int memory 
        Int disk_space 
        Int num_threads 
        }
    call plink2AF {
        input:
            pvar = pvar, 
            psam = psam,
            pgen = pgen,
            prefix = prefix,
            sample_list = sample_list,
            memory = memory, 
            disk_space = disk_space, 
            num_threads = num_threads
        }
    output {
        File PlinkAF = plink2AF.PlinkAF
        }
    }

    task plink2AF {
        input {
            File pvar 
            File pgen 
            File psam
            File sample_list
            String prefix
            Int memory 
            Int disk_space 
            Int num_threads
        }

        command <<<
        first=$(head -n1 ~{sample_list})
        if [[ "$first" =~ [a-zA-Z] ]]; then
          tail -n +2 ~{sample_list} | awk '{print $1}' > keep_plink.txt
        else
          awk '{print $1}' ~{sample_list} > keep_plink.txt
        fi

        plink2 \
            --pvar "~{pvar}" \
            --psam "~{psam}" \
            --pgen "~{pgen}" \
            --keep keep_plink.txt \
            --freq --out "~{prefix}" 
        >>> 
        
        runtime {
            docker: "quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0"
            memory: "~{memory}GB"
            disks: "local-disk ~{disk_space} HDD"
            cpu: "~{num_threads}"

        }

        output {
            File PlinkAF = "~{prefix}.afreq" 
        }
    }
