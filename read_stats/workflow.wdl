version 1.0

workflow read_stats {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Simple stats about reads in one gzipped FASTQ file. Simplified from HPP workflow https://github.com/human-pangenomics/hpp_production_workflows/blob/master/QC/wdl/tasks/read_stats.wdl (one file input, must be gzipped FASTQ, but all in one task)"
    }
    input {
        File READS
        Int histogramMinLength=0
        Int histogramMaxLength=0
    }

    call readStats {
        input:
        readsFile=READS,
        histogramMinLength=histogramMinLength,
        histogramMaxLength=histogramMaxLength
    }
    
    output {
        File output_tar = readStats.outputTarball
        File output_report = readStats.reportTsv
    }
}

task readStats {
    input {
        File readsFile
        Int histogramMinLength = 0
        Int histogramMaxLength = 0
        Int memSizeGB = 2
        Int threadCount = 2
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

        ln -s ~{readsFile}
                
        FILE=$(basename ~{readsFile})
        FAI_OUTPUT="$FILE.fai"

        zcat $FILE | awk '{if(NR%4==2) print length($1)}' > $FAI_OUTPUT
        
        # get output name
        OUTPUT=$(basename ~{readsFile} | sed -E 's/.(fastq.gz|fq.gz)*$//')

        # hist parameters
        if [[ ~{histogramMinLength} -eq 0 && ~{histogramMaxLength} -eq 0 ]] ; then
            HIST_PARAM="--hist_auto_bounds"
        else
            HIST_PARAM="--hist_min ~{histogramMinLength} --hist_max ~{histogramMaxLength}"
        fi

        # sketch
        fai_read_stats.py -i $FAI_OUTPUT -o $OUTPUT $HIST_PARAM

        # for consolidation
        cat $OUTPUT/report.tsv | sed "s/^/$OUTPUT\t/" > $OUTPUT.report.tsv

        # rename and tar
        ls $OUTPUT/ | xargs -n 1 -I{} mv $OUTPUT/{} $OUTPUT/$OUTPUT.{}
        tar czvf $OUTPUT.tar.gz $OUTPUT/

	>>>

	output {
		File outputTarball = glob("*.tar.gz")[0]
		File reportTsv = glob("*.report.tsv")[0]
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "humanpangenomics/fai_read_stats:latest"
        preemptible: 1
    }
}
