version 1.0

workflow add_unmapped_gam_bam {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Add missing alignment in a BAM from a GAM"
    }
    input {
        File GAM
        File BAM
        File SCRIPT
    }

    call extractSortedReadNameGam {
        input:
        gam=GAM
    }
    
    call extractSortedReadNameBam {
        input:
        bam=BAM
    }

    call listMissingReadnames {
        input:
        all_reads=extractSortedReadNameGam.readNames,
        aligned_reads=extractSortedReadNameBam.readNames
    }

    call subsetIntoUmappedBam {
         input:
         gam=GAM,
         read_names=listMissingReadnames.readNames,
         script=SCRIPT
    }

    call mergeBams {
        input:
        bam1=BAM,
        bam2=subsetIntoUmappedBam.bam
    }
    
    output {
        File bam = mergeBams.merged_bam
        File bam_index = mergeBams.merged_bam_index
    }
}

task extractSortedReadNameGam {
    input {
        File gam
        Int memSizeGB = 4
        Int threadCount = 3
        Int diskSizeGB = 2*round(size(gam, "GB")) + 10
    }

	command <<<
 set -eux -o pipefail
 touch empty.gam
 vg gamcompare -T ~{gam} empty.gam | sed 1d | cut -f 4 | sort -S 1G | gzip > gam.rnames.txt.gz
	>>>

	output {
		File readNames = "gam.rnames.txt.gz"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/vgteam/vg:v1.44.0"
        preemptible: 1
    }
}


task extractSortedReadNameBam {
    input {
        File bam
        Int memSizeGB = 4
        Int threadCount = 3
        Int diskSizeGB = 2*round(size(bam, "GB")) + 10
    }

	command <<<
        set -eux -o pipefail

        samtools view ~{bam} | cut -f 1 | sort -S 1G | gzip > bam.rnames.txt.gz
	>>>

	output {
		File readNames = "bam.rnames.txt.gz"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/vgteam/vg:v1.44.0"
        preemptible: 1
    }
}

task listMissingReadnames {
    input {
        File all_reads
        File aligned_reads
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 50*round(size(all_reads, "GB")) + 10
    }

	command <<<
        set -eux -o pipefail
        
        gunzip -c ~{all_reads} > all.txt
        gunzip -c ~{aligned_reads} > aligned.txt
        comm -23 all.txt aligned.txt | gzip > missing.rnames.txt.gz
	>>>

	output {
		File readNames = "missing.rnames.txt.gz"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/vgteam/vg:v1.44.0"
        preemptible: 1
    }
}

task subsetIntoUmappedBam {
    input {
        File gam
        File read_names
        File script
        Int memSizeGB = 4
        Int threadCount = 3
        Int diskSizeGB = 5*round(size(gam, "GB") + size(read_names, "GB")) + 20
    }

	command <<<
        set -eux -o pipefail
        
        gunzip -c ~{read_names} > readnames.txt

        vg view -a ~{gam} | jq -r "[.name,.sequence] | @tsv" | python3 ~{script} -r readnames.txt | samtools view -bh - > unmapped.bam
	>>>

	output {
		File bam = "unmapped.bam"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/vgteam/vg:v1.44.0"
        preemptible: 1
    }
}

task mergeBams {
    input {
        File bam1
        File bam2
        Int memSizeGB = 4
        Int threadCount = 2
        Int diskSizeGB = 5*round(size(bam1, "GB") + size(bam2, "GB")) + 20
    }

	command <<<
        set -eux -o pipefail
        samtools merge -O BAM merged.bam ~{bam1} ~{bam2}

        samtools index -b merged.bam merged.bam.bai
	>>>

	output {
		File merged_bam = "merged.bam"
		File merged_bam_index = "merged.bam.bai"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
        preemptible: 1
    }
}
