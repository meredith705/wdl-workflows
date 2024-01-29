version 1.0

workflow chimera_stats {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Estimate the amount of chimeric (long) reads. Uses https://github.com/rlorigro/Liger2LiGer."
    }
    input {
        File? FASTQ_FILE                  ## input reads in gzipped FASTQ
        File? REFERENCE_FILE              ## reference fasta
        File? BAM_FILE                    ## input reads in gzipped FASTQ
    }

    if(defined(FASTQ_FILE) && defined(REFERENCE_FILE)){
        call runLiger2LiGerFromFastq {
            input:
            readsFile=FASTQ_FILE,
            referenceFile=REFERENCE_FILE
        }
    }

    if(defined(BAM_FILE)){
        call runLiger2LiGerFromBam {
            input:
            bamFile=BAM_FILE
        }
    }

    File outTar = select_first([runLiger2LiGerFromFastq.outputTarball, runLiger2LiGerFromBam.outputTarball])
    File outReport = select_first([runLiger2LiGerFromFastq.reportCsv, runLiger2LiGerFromBam.reportCsv])
    
    output {
        File output_tar = outTar
        File output_report = outReport
    }
}

task runLiger2LiGerFromFastq {
    input {
        File? readsFile
        File? referenceFile
        Int memSizeGB = 50
        Int threadCount = 16
    }

    Int diskSizeGB = round(size(readsFile, 'G') + size(referenceFile, 'G')) * 5 + 20
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

        python3 /build/Liger2LiGer/scripts/evaluate_chimeras.py --ref ~{referenceFile} --fastq ~{readsFile} --n_threads ~{threadCount}

        OUTPREF=$(basename ~{readsFile} | sed -E 's/.(fastq.gz|fq.gz)*$//')

        ## rename report
        mv results_*csv ${OUTPREF}.liger2liger_results.csv

        ## tar some outputs
        mkdir -p ${OUTPREF}.liger2liger_output
        mv */*.png ${OUTPREF}.liger2liger_output/
        mv */*.txt ${OUTPREF}.liger2liger_output/
        tar -czvf ${OUTPREF}.liger2liger_output.tar.gz ${OUTPREF}.liger2liger_output/

	>>>

	output {
		File outputTarball = glob("*liger2liger_output.tar.gz")[0]
		File reportCsv = glob("*liger2liger_results.csv")[0]
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/liger2liger@sha256:ec76cec82c11323806b4e4681cbce7894a334c6338e664d1e49eaefacd757731"
        preemptible: 1
    }
}

task runLiger2LiGerFromBam {
    input {
        File? bamFile
        Int memSizeGB = 4
        Int threadCount = 1
    }

    String outprefix = basename(select_first([bamFile]), '.bam')
    Int diskSizeGB = 2*round(size(bamFile, 'G')) + 20
    
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

        ln -s ~{bamFile} reads.bam
        
        filter_chimeras_from_alignment -i reads.bam
        
        python3 /build/Liger2LiGer/scripts/generate_chimer_stats.py \
                -i reads.chimer_lengths.txt,reads.non_chimer_lengths.txt \
                -o liger2liger_output
        
        ## rename report
        mv liger2liger_output/results_*csv ~{outprefix}.liger2liger_results.csv
        
        ## tar some outputs
        mkdir -p ~{outprefix}.liger2liger_output
        mv liger2liger_output/*.png ~{outprefix}.liger2liger_output/
        mv *.txt ~{outprefix}.liger2liger_output/
        tar czvf ~{outprefix}.liger2liger_output.tar.gz ~{outprefix}.liger2liger_output/
	>>>

	output {
		File outputTarball = "~{outprefix}.liger2liger_output.tar.gz"
		File reportCsv = "~{outprefix}.liger2liger_results.csv"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/liger2liger@sha256:ec76cec82c11323806b4e4681cbce7894a334c6338e664d1e49eaefacd757731"
        preemptible: 1
    }
}
