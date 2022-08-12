version 1.0

workflow align_stats {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Compute the distribution of alignment identity using wambam: https://github.com/rlorigro/wambam."
    }

    parameter_meta {
        BAM_FILE: "BAM file from running minimap2 with the --eqx flag. Either provide this file or both a FASTQ_FILE and REFERENCE_FILE."
        FASTQ_FILE: "Reads in a gzipped FASTQ file. Either provide this file and REFERENCE_FILE, or a BAM_FILE."
        REFERENCE_FILE: "FASTA file for the reference genome. Can be gzipped. Either provide this file and FASTQ_FILE, or a BAM_FILE."
    }

    input {
        File? BAM_FILE
        File? FASTQ_FILE
        File? REFERENCE_FILE
    }

    if(!defined(BAM_FILE) && defined(FASTQ_FILE) && defined(REFERENCE_FILE)){
        call runMinimap2 {
            input:
            readsFile=FASTQ_FILE,
            referenceFile=REFERENCE_FILE
        }
    }

    File cur_bam_file = select_first([BAM_FILE, runMinimap2.bam])
    
    call runWambam {
        input: bamFile=cur_bam_file
    }
    
    output {
        File identity_dist_csv = runWambam.outputIdDist
        File? bam = runMinimap2.bam
        File? bam_index = runMinimap2.bam_index
    }
}

task runWambam {
    input {
        File bamFile
        Int memSizeGB = 10
    }

    Int diskSizeGB = round(2*size(bamFile, "GB")) + 50

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
        docker: "quay.io/jmonlong/wambam@sha256:55bd4b233e6cb104a70919ae38838ac5914dfb7bb5e7c5ef94b39eb5c698bbd0"
        preemptible: 1
    }
}

task runMinimap2 {
    input {
        File? readsFile
        File? referenceFile
        String preset = "map-ont"
        Int kSize = 15
	    Boolean useMd = true
        Int memSizeGB = 128
        Int threadCount = 64
    }

    Int diskSizeGB = 10 * round(size(readsFile, "GB") + size(referenceFile, "GB")) + 50
    Int threadMinimap = if threadCount < 5 then threadCount else threadCount - 4
    Int threadSort = if threadCount < 5 then 1 else 4
    Int memSort = if memSizeGB < 8 then 1 else 4
    
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

        minimap2 -x ~{preset} -K 5G ~{true="--MD" false="" useMd} -a -c --eqx -t ~{threadMinimap} -k ~{kSize} ~{referenceFile} ~{readsFile} | samtools sort -@ ~{threadSort} -m ~{memSort}G > $OUTPREF.bam
        samtools index -@ ~{threadCount} $OUTPREF.bam $OUTPREF.bam.bai
	>>>

	output {
		File bam = glob("*.bam")[0]
		File bam_index = glob("*.bam.bai")[0]
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "mkolmogo/card_mapping"
        preemptible: 1
    }
}
