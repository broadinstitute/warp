version 1.0

workflow HLAGenotyping {

    String pipeline_version = "aou_9.0.0"

    input {
        String gatk_docker
        String hlahd_docker
        String polysolver_docker
        String optitype_docker

        File original_bam       # this can be a BAM or CRAM
        File original_bam_idx
        File ref_fasta
        File ref_fai
        File ref_dict
        File hla_intervals

        File convert_alleles_python_script
        File count_two_field_alleles_python_script
        File hla_groups_file

        String? gcs_project_for_requester_pays
        # WDL version 1.0 does not have an empty Optional literal
        # such a literal is very useful because Terra has a bug where whenever a data table is updated, empty values
        # silently and invisibly get converted to empty strings "".  Thus it is useful to recognize empty strings and
        # declare empty Optionals.  The only way to do this in WDL 1.0 is to get an empty Optional as a variable from the
        # workflow inputs.  These inputs should NEVER be filled in!!!!!
        File? EMPTY_STRING_HACK
    }


    call MakeHLAOnlyBamsAndFastqs {
        input:
            gatk_docker = gatk_docker,
            gcs_project_for_requester_pays = if select_first([gcs_project_for_requester_pays, ""]) == "" then EMPTY_STRING_HACK else gcs_project_for_requester_pays,
            original_bam = original_bam,
            original_bam_idx = original_bam_idx,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            hla_intervals = hla_intervals
    }

    call HLAHD {
        input:
            hlahd_docker = hlahd_docker,
            first_end_fastq = MakeHLAOnlyBamsAndFastqs.first_end_fastq,
            second_end_fastq = MakeHLAOnlyBamsAndFastqs.second_end_fastq,
            convert_alleles_python_script = convert_alleles_python_script,
            count_two_field_alleles_python_script = count_two_field_alleles_python_script,
            hla_groups_file = hla_groups_file
    }



    if (HLAHD.two_field_count > 0) {
        call Polysolver {
            input:
                polysolver_docker = polysolver_docker,
                hla_bam = MakeHLAOnlyBamsAndFastqs.hla_bam,
                hla_bam_idx = MakeHLAOnlyBamsAndFastqs.hla_bam_idx
        }

        call Optitype {
            input:
                optitype_docker = optitype_docker,
                first_end_fastq = MakeHLAOnlyBamsAndFastqs.first_end_fastq,
                second_end_fastq = MakeHLAOnlyBamsAndFastqs.second_end_fastq
        }

        call Consensus {
	        input:
		        hla_hd_result = HLAHD.converted_result,
                polysolver_result = Polysolver.result,
                optitype_result = Optitype.result
        }
    }

    output {
        File hlahd_raw_result = HLAHD.raw_result
        File hlahd_converted_result = HLAHD.converted_result
        Int hlahd_two_field_count = HLAHD.two_field_count
        Int hlahd_overriden = select_first([Consensus.hlahd_overridden, 0])
        File consensus = select_first([Consensus.consensus, HLAHD.converted_result])

        File? optitype_result = Optitype.result
        File? polysolver_result = Polysolver.result
    }
}

task MakeHLAOnlyBamsAndFastqs {
    input {
        String gatk_docker
        String? gcs_project_for_requester_pays
        File original_bam       # this can be a BAM or CRAM
        File original_bam_idx
        File ref_fasta          # GATK PrintReads requires a reference for CRAMs
        File ref_fai
        File ref_dict
        File hla_intervals

        Int cpu = 2
        Int num_threads = 4
        Int mem_gb = 4
        Int disk_gb = 100
        Int boot_disk_gb = 10
        Int max_retries = 0
        Int preemptible = 1
    }

    parameter_meta{
        hla_intervals: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fai: {localization_optional: true}
        ref_dict: {localization_optional: true}
        original_bam: {localization_optional: true}
        original_bam_idx: {localization_optional: true}
    }

    command <<<
        # this command also produces the accompanying index hla.bai
        # the PairedReadFilter is necessary for SamtoFastq to succeed
        gatk PrintReads -R ~{ref_fasta} -I ~{original_bam} -L ~{hla_intervals} -O hla-unsorted.bam \
            ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}


        echo "We are running ValidateSamFile on the output of PrintReads:"
        gatk ValidateSamFile -I hla-unsorted.bam

        samtools sort -@ ~{num_threads} hla-unsorted.bam > hla.bam

        # using gatk instead of samtools for indexing avoids ERROR:INVALID_INDEX_FILE_POINTER in the output
        # of ValidateSamFile.  I'm not sure what that means or if it's important, but might as well not have it.
        gatk BuildBamIndex -I hla.bam
        #samtools index -@ ~{num_threads} hla.bam > hla.bai

        echo "We are running ValidateSamFile on the output of samtools sorting and re-indexing"
        gatk ValidateSamFile -I hla.bam

        # The "*" MUST be in quotes.  -T "*" indicates that all tags are copied to output.
        samtools fastq -@ ~{num_threads} -n -T "*" -0 /dev/null -1 first_end.fq -2 second_end.fq hla.bam
    >>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: boot_disk_gb
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    output {
        File hla_bam = "hla.bam"
        File hla_bam_idx = "hla.bai"
        File first_end_fastq = "first_end.fq"
        File second_end_fastq = "second_end.fq"
    }
}

task HLAHD {
    input {
        String hlahd_docker
        File first_end_fastq
        File second_end_fastq

        File convert_alleles_python_script
        File count_two_field_alleles_python_script
        File hla_groups_file

        Int min_read_length = 75
        Float trim_ratio = 0.95

        Int cpu = 2
        Int num_threads = 4

        Int mem_gb = 20
        Int disk_gb = 100
        Int boot_disk_gb = 25
        Int max_retries = 0
        Int preemptible = 0
    }

    command <<<
        mkdir OUTPUT

        hlahd.sh -t ~{num_threads} -m ~{min_read_length} -c ~{trim_ratio} \
            -f /HLA/hlahd.1.7.1/freq_data/ \
            ~{first_end_fastq} ~{second_end_fastq} \
            /HLA/hlahd.1.7.1/HLA_gene.split.3.32.0.txt \
            /HLA/hlahd.1.7.1/dictionary SAMPLE_ID \
            OUTPUT

        mv OUTPUT/SAMPLE_ID/result/SAMPLE_ID_final.result.txt ./raw1.txt

        # HLA-HD alleles are in the form HLA-A*23:01:01
        # here we remove the "HLA-" from every allele, remove lines that say Couldn't read result file,
        # and convert "Not typed" (with a pesky space) to "NA"
        sed 's/HLA-//g' raw1.txt | grep -v "Couldn't read result" | sed 's/Not typed/NA/g' > raw2.txt

        # HLA-HD emits a '-' in the second allele to denote a homozygous genotype.  We expand this out.
        # also, it sometimes emits a third and fourth allele, which are weaker guesses that we discard
        # finally, we sort in lexicographical order
        touch raw3.txt
        while read gene allele1 allele2 other_alleles; do
            printf "${gene}\t" >> raw3.txt
            if [[ ${allele2} == "-" ]]; then
               printf "${allele1}\t${allele1}\n" >> raw3.txt   # copy allele1 for homozygous
            elif [[ ${allele1} < ${allele2} ]]; then
               printf "${allele1}\t${allele2}\n" >> raw3.txt
            else
               printf "${allele2}\t${allele1}\n" >> raw3.txt
            fi
         done < raw2.txt

        python3 ~{convert_alleles_python_script} raw3.txt ~{hla_groups_file} > hlahd_converted.txt

        python3 ~{count_two_field_alleles_python_script} hlahd_converted.txt > two_field_count.txt
    >>>

    runtime {
        docker: hlahd_docker
        bootDiskSizeGb: boot_disk_gb
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    output {
        File raw_result = "raw3.txt"
        File converted_result = "hlahd_converted.txt"
        Int two_field_count = read_int("two_field_count.txt")
    }
}

task Polysolver {
    input {
        String polysolver_docker
        String polysolver_reference = "hg38"
        File hla_bam
        File hla_bam_idx

        Int cpu = 2
        Int num_threads = 4
        Int mem_gb = 20
        Int disk_gb = 100
        Int boot_disk_gb = 25
        Int max_retries = 0
        Int preemptible = 0
    }

    command <<<

        # polysolver (positional) args are, in order:
        #   -bam: path to the BAM file to be used for HLA typing
        #	-race: ethnicity of the individual (Caucasian, Black, Asian or Unknown)
        #	-includeFreq: flag indicating whether population-level allele frequencies should be used as priors (0 or 1)
        #	-build: reference genome used in the BAM file (hg18, hg19 or hg38)
        #	-format: fastq format (STDFQ, ILMFQ, ILM1.8 or SLXFQ; see Novoalign documentation)
        #	-insertCalc: flag indicating whether empirical insert size distribution should be used in the model (0 or 1)
        #	-outDir: output directory

        mkdir output

        /home/polysolver/scripts/shell_call_hla_type ~{hla_bam} Unknown 1 ~{polysolver_reference} STDFQ 0 output

        echo "Here is ls of the output directory:"
        ls output

        cd /cromwell_root
        mv output/* .

        mv winners.hla.txt raw1.txt

        echo "RAW1"
        cat raw1.txt

        # Polysover output looks like
        # HLA-A	hla_a_23_01_01	hla_a_34_02_01
        # HLA-B	hla_b_44_03_01	hla_b_45_01
        # HLA-C	hla_c_04_01_01_04	hla_c_06_02_01_01
        # we convert HLA-A to just A, delete "hla_", convert a_23_01_01 to A*23:01:01, and remove the fourth field

        sed 's/HLA-//g' raw1.txt | sed 's/hla_//g' | tr '[:lower:]' '[:upper:]' > raw2.txt
        sed 's/\([A-C]\)_/\1\*/g' raw2.txt | sed 's/_/:/g' > raw3.txt
        sed 's/\([0-9]*:[0-9]*:[0-9]*\):[0-9]*/\1/g' raw3.txt > raw4.txt

        echo "RAW2"
        cat raw2.txt

        echo "RAW3"
        cat raw3.txt

        echo "RAW4"
        cat raw4.txt

        # finally, we sort in lexicographical order
        touch sorted.txt
        while read gene allele1 allele2; do
            printf "${gene}\t" >> sorted.txt
            if [[ ${allele1} < ${allele2} ]]; then
               printf "${allele1}\t${allele2}\n" >> sorted.txt
            else
               printf "${allele2}\t${allele1}\n" >> sorted.txt
            fi
         done < raw4.txt


    >>>

    runtime {
        docker: polysolver_docker
        bootDiskSizeGb: boot_disk_gb
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    output {
        File result = "sorted.txt"
        File check = "check.status.out.txt"
        File counts1 = "counts1.R0k6"
        File counts2 = "counts2.R0k6"
    }
}

task Optitype {
    input {
        String optitype_docker
        File first_end_fastq
        File second_end_fastq

        Int mem_gb = 20
        Int disk_gb = 100
        Int boot_disk_gb = 25
        Int max_retries = 0
        Int preemptible = 0
        Int cpu = 1
    }

    command <<<
        python /usr/local/bin/OptiType/OptiTypePipeline.py --input ~{first_end_fastq} ~{second_end_fastq} --dna --verbose --outdir ./output

        ls ./output

        #Optitype puts output not in the specificed outdir dir directly but rather in a time-stamped subfolder
        echo "Here are tsvs"
        find ./output -type f -name "*.tsv"

        echo "here is everything"
        find ./output -name "*"

        find ./output -type f -name "*.tsv" | while read file; do cp $file .; done

        echo "Here is all tsv"
        ls *.tsv

        mv *.tsv optitype_raw.tsv

        # Optitype format is two lines of output
        # A1	A2	B1	B2	C1	C2	Reads	Objective
        # 0	A*23:01	A*34:02	B*45:04	B*44:03	C*06:02	C*04:01	845.0	806.9750000000003

        # we convert it to
        # A A*23:01	A*34:02
        # B B*45:04	B*44:03
        # C C*06:02	C*04:01

        tail -n +2 optitype_raw.tsv | cut -f 2-7 | \
            while read a1 a2 b1 b2 c1 c2; do printf "A\t${a1}\t${a2}\nB\t${b1}\t${b2}\nC\t${c1}\t${c2}\n"; done > optitype.tsv

        # finally, we sort in lexicographical order
        touch sorted.txt
        while read gene allele1 allele2; do
            printf "${gene}\t" >> sorted.txt
            if [[ ${allele1} < ${allele2} ]]; then
               printf "${allele1}\t${allele2}\n" >> sorted.txt
            else
               printf "${allele2}\t${allele1}\n" >> sorted.txt
            fi
         done < optitype.tsv
    >>>

    runtime {
        docker: optitype_docker
        bootDiskSizeGb: boot_disk_gb
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }

    output {
        File result = "sorted.txt"
    }
}

task Consensus {
	input {
		File hla_hd_result
        File polysolver_result
        File optitype_result
		Int disk_space = 10
	}

	command <<<
        # from previous tasks, all the files are in a standard format of eg:
        # A   A*15:01:02        A*17:02:01
        # B   B*3:01:02         B*3:01:02   # (note that for homozygous we have converted dashes in HLA-HD)
        # etc, with Polysolver and Optitype only ioncluding A, B, C and HLA-HD going further.
        # Optitype has two-field resolution, Polysolver has 4, HLA-HD usually has 3 but sometimes 2
        # the two alleles for each gene are sorted lexicographically

        head -n 3 ~{hla_hd_result} > HLAHD_ABC
        tail -n +4 ~{hla_hd_result} > HLA_HD_REST

        # make two-field versions of everything
        for file in HLAHD_ABC ~{polysolver_result} ~{optitype_result}; do
            sed 's/\([0-9]*:[0-9]*\):[0-9]*/\1/g' ${file} > ${file}_TWO_FIELD

            echo "original file ${file}:"
            cat ${file}

            echo "two-field file ${file}_TWO_FIELD:"
            cat ${file}_TWO_FIELD
        done

        # paste the everything together
        paste ~{polysolver_result}_TWO_FIELD ~{polysolver_result} HLAHD_ABC_TWO_FIELD HLAHD_ABC ~{optitype_result}_TWO_FIELD > PASTED

        echo "Pasted:"
        cat PASTED

        hlahd_overridden=0
        touch consensusABC
        while read gene1 ps1 ps2 gene2 ps3f1 ps3f2 gene3 hlahd1 hlahd2 gene4 hlahd3f1 hlahd3f2 gene5 opt1 opt2; do
             ps=${ps1}${ps2}
             hlahd=${hlahd1}${hlahd2}
             opt=${opt1}${opt2}

            if [[ ${hlahd} != ${ps} && ${ps} == ${opt} ]]; then
                ((hlahd_overridden++))
                printf "${gene1}\t${ps3f1}\t${ps3f2}\n" >> consensusABC
            else
                printf "${gene1}\t${hlahd3f1}\t${hlahd3f2}\n" >> consensusABC
            fi
      done < PASTED

      echo ${hlahd_overridden} > hlahd_overridden.txt
      cat consensusABC HLA_HD_REST > consensus.txt
	>>>

    runtime {
        docker: "continuumio/anaconda:latest"
        disks: "local-disk " + disk_space + " SSD"
    }

    output {
        File consensus = "consensus.txt"
        Int hlahd_overridden = read_int("hlahd_overridden.txt")
    }
}
