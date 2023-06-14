version 1.0

workflow vcf_stats {

    meta {
        author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Extract some stats from a VCF: phase block N50, number of het/hom SVs"
    }

    parameter_meta {
        VCF: "Input VCF. Can be gzipped/bgzipped."
    }
    
    input {
        File VCF
    }

    call run_vcf_stats {
        input:
        vcf=VCF
    }
        
    output {
        File vcf_stats = run_vcf_stats.tsv
    }
}

task run_vcf_stats {
    input {
        File vcf
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 2*round(size(vcf, "GB")) + 10
    }

    String basen = sub(sub(basename(vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
    command <<<
    set -eux -o pipefail

    python3 /opt/scripts/get_vcf_stats.py -v ~{vcf} -o ~{basen}.vcf_stats.tsv
    >>>
    
    output {
        File tsv = "~{basen}.vcf_stats.tsv"
    }
    
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/vcf_stats"
        preemptible: 1
    }
}

