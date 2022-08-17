version 1.0

workflow read_subset_qscore {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Subset a fastq for reads with a minimum basecalling qscore."
    }
    parameter_meta {
        FASTQ_FILE: "Reads in a gzipped FASTQ file. "
        SUMMARY_FILE: "Guppy summary (gzipped TSV file). "
        MIN_QSCORE: "keep only reads with a qscore greater or equal to thie value. Default: 10"
    }

    input {
        File FASTQ_FILE
        File SUMMARY_FILE
        Int MIN_QSCORE=10
    }

    call subsetReads {
        input:
        fastq_file=FASTQ_FILE,
        summary_file=SUMMARY_FILE,
        min_qscore=MIN_QSCORE
    }
    
    output {
        File output_fastq = subsetReads.fastq
    }
}

task subsetReads {
    input {
        File fastq_file
        File summary_file
        Int min_qscore = 10
        Int memSizeGB = 4
        Int threads = 4
    }

    Int diskSizeGB = 3 * round(size(fastq_file, 'G')) + 10
    Int threadsGzip = if threads > 1 then threads - 1 else 1
	command <<<
        set -eux -o pipefail

        python3 /scripts/subset_reads_by_qscore.py -f ~{fastq_file} -s ~{summary_file} -q ~{min_qscore} | pigz -p {threadsGzip} > reads_min~{min_qscore}.fastq.gz
	>>>

	output {
		File fastq = "reads_min~{min_qscore}.fastq.gz"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/read_subset_qscore:latest"
        preemptible: 1
    }
}
