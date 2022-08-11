version 1.0

workflow align_stats {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Compute the distribution of alignment identity using wambam: https://github.com/rlorigro/wambam"
    }
    input {
        File READS                       ## input reads in gzipped FASTQ
        File REFERENCE_FILE              ## reference fasta
    }

    call runMinimap2 {
        input:
        readsFile=READS,
        referenceFile=REFERENCE_FILE
    }

    call runWambam {
        input:
        bamFile=runMinimap2.bam
    }
    
    output {
        File output_identity_dist = runWambam.outputIdDist
    }
}

task runWambam {
    input {
        File bamFile
        Int memSizeGB = 10
        Int diskSizeGB = 200
    }

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

        wam -i ~{bamFile} -o wambam_results
	>>>

	output {
		File outputIdDist = "wambam_results/identity_distribution.csv"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: 1
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/wambam@sha256:4163e0bd0f7969a960d5bc777ab0fe67e9cecf423c7e9cf9781aded37889d9f6"
        preemptible: 1
    }
}

task runMinimap2 {
    input {
        File readsFile
        File referenceFile
        String preset = "map-ont"
        Int kSize = 15
        Int memSizeGB = 40
        Int threadCount = 16
        Int diskSizeGB = 500
    }

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

        OUTPREF=$(basename ~{readsFile} | sed -E 's/.(fastq.gz|fq.gz)*$//')

        minimap2 -x ~{preset} -a -c --eqx -t ~{threadCount} -k ~{kSize} ~{referenceFile} ~{readsFile} | samtools view -hB > $OUTPREF.bam
	>>>

	output {
		File bam = glob("*.bam")[0]
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "mkolmogo/card_mapping"
        preemptible: 1
    }
}
