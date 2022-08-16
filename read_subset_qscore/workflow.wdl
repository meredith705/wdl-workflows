version 1.0

workflow read_subset_qscore {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Subset a fastq for reads with a minimum basecalling qscore. Both fastq and summary inputs should have the reads in the same order. They should if were the output of the same Guppy run. "
    }
    parameter_meta {
        FASTQ_FILE: "Reads in a gzipped FASTQ file. Should be in the same order as in the SUMMARY_FILE"
        SUMMARY_FILE: "Guppy summary (gzipped TSV file). Should be in the same order as in the FASTQ_FILE"
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
        Int memSizeGB = 2
    }

    Int diskSizeGB = 3 * round(size(fastq_file, 'G')) + 10
	command <<<
        set -eux -o pipefail

        python3 /scripts/subset_reads_by_qscore.py -f ~{fastq_file} -s ~{summary_file} -q ~{min_qscore} -o reads_min~{min_qscore}.fastq.gz
	>>>

	output {
		File fastq = "reads_min~{min_qscore}.fastq.gz"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: 1
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "jmonlong/read_subset_qscore:latest"
        preemptible: 1
    }
}
