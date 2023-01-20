version 1.0

workflow indexBam {
    meta {
	author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Index a BAM file" 
    }
    input {
        File BAM_FILE
    }

    call bamIndex {
        input:
        bam=BAM_FILE
    }
    
    output {
        File bam_index = bamIndex.bam_index
    }
}

task bamIndex {
    input {
        File bam
        Int disk = 2 * round(size(bam, 'G')) + 20
    }
    String outpref = basename(bam)
    command <<<
    set -eux -o pipefail
    
    samtools index -b ~{bam} ~{outpref}.bai
    >>>
    output {
        File bam_index= "~{outpref}.bai"
    }
    runtime {
        docker: "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk " + disk + " SSD"
        preemptible: 1
    }
}

