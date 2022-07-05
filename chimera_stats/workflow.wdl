version 1.0

workflow chimera_stats {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Estimate the amount of chimeric (long) reads. Uses https://github.com/rlorigro/Liger2LiGer."
    }
    input {
        File READS                       ## input reads in gzipped FASTQ
        File REFERENCE_FILE              ## reference fasta
    }

    call runLiger2LiGer {
        input:
        readsFile=READS,
        referenceFile=REFERENCE_FILE
    }
    
    output {
        File output_tar = runLiger2LiGer.outputTarball
        File output_report = runLiger2LiGer.reportCsv
    }
}

task runLiger2LiGer {
    input {
        File readsFile
        File referenceFile
        Int memSizeGB = 10
        Int threadCount = 16
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

        python3 /build/Liger2LiGer/scripts/evaluate_chimeras.py --ref ~{referenceFile} --fastq ~{readsFile} --n_threads ~{threadCount}

        OUTPREF=$(basename ~{readsFile} | sed -E 's/.(fastq.gz|fq.gz)*$//')

        ## rename report
        mv results_*csv ${OUTPREF}.liger2liger_results.csv

        ## tar some outputs
        mkdir -p ${OUTPREF}.liger2liger_output
        mv */*chimer_distribution_coverage.png ${OUTPREF}.liger2liger_output/
        mv */*chimer_distribution_percent.png ${OUTPREF}.liger2liger_output/
        mv */*chimer_distribution_raw.png ${OUTPREF}.liger2liger_output/
        mv */*.chimeric_reads.txt ${OUTPREF}.liger2liger_output/
        tar czvf ${OUTPREF}.liger2liger_output.tar.gz ${OUTPREF}.liger2liger_output/

	>>>

	output {
		File outputTarball = glob("*liger2liger_output.tar.gz")[0]
		File reportCsv = glob("*liger2liger_results.csv")[0]
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/liger2liger@sha256:a6b8b50482d5dd6af2e2d6ab4ba0ae78c661b5339d5f26e9ff224e10f3d9abdd"
        preemptible: 1
    }
}
