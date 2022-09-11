version 1.0

workflow quast {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Run QUAST to get QC metrics on an assembly slightly tweaked from https://dockstore.org/workflows/github.com/human-pangenomics/hpp_production_workflows/Quast:master?tab=info to allow for an input array of fasta files (that will be concatenated)."
    }

    parameter_meta {
        AMB_FA: "Assembly file(s). One or multiple fasta files."
        REF_FA: "Reference genome fasts (optional)"
    }

    input {
        Array[File] AMB_FA
        File? REF_FA
    }

    call quast {
        input:
        assemblyFasta=AMB_FA,
        referenceFasta=REF_FA
    }
    
    output {
		File quastTarball = quast.outputTarball
		File quastSummary = quast.outputSummary
    }
}


task quast {
    input {
        Array[File] assemblyFasta
        File? referenceFasta
        String? label = "assembly"
        String extraArguments="--large"
        Int memSizeGB = 64
        Int threadCount = 16
        String dockerImage = "tpesout/hpp_quast:latest"
    }

    Int diskSizeGB = round(5*(size(assemblyFasta, "GB") + size(referenceFasta, "GB"))) + 50

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # initilization
        rm -f ~{label}.fasta
        for $ASM_FULL_PATH in ~{sep=" " assemblyFasta}
        do
            ASM_FILENAME=$(basename -- "$ASM_FULL_PATH")
            if [[ $ASM_FILENAME =~ \.gz$ ]]; then
                gunzip -c $ASM_FULL_PATH >> ~{label}.fasta
            else
                cat $ASM_FULL_PATH >> ~{label}.fasta
            fi
        done

        # init quast command
        cmd=(python /opt/quast/quast-5.0.2/quast-lg.py )
        cmd+=( -t ~{threadCount} )
        cmd+=( -o ~{label}.quast )

        # include reference fasta if supplied
        if [[ -f "~{referenceFasta}" ]]; then
            REF_FILENAME=$(basename -- "~{referenceFasta}")
            if [[ $REF_FILENAME =~ \.gz$ ]]; then
                cp ~{referenceFasta} .
                gunzip $REF_FILENAME
                REF_FILENAME="${REF_FILENAME%.gz}"
            else
                ln -s ~{referenceFasta}
            fi
            cmd+=( -r $REF_FILENAME )
        fi

        # include extra arguments if supplied
        if [[ ! -z "~{extraArguments}" ]] ; then
            cmd+=( ~{extraArguments} )
        fi

        # finalize command
        cmd+=( ~{label}.fasta )

        # run command
        "${cmd[@]}"

        # save output
        tar czvf ~{label}.quast.tar.gz ~{label}.quast

	>>>
	output {
		File outputTarball = "~{label}.quast.tar.gz"
		File outputSummary = "~{label}.quast/report.txt"
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
