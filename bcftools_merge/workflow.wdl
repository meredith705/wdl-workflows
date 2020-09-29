version 1.0

workflow bcftools_merge {
    meta {
	author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Merge VCFs using bcftools. Multi-allelic records are split and the merged VCF is normalized (requires a reference fasta)."
    }
    input {
        Array[File] VCF_LIST
        File REF_FASTA
        File REF_FASTA_IDX
        Int CORES = 1
        Int DISK = 100
        Int MEM = 10
        Int PREEMPTIBLE = 0
    }

    call merge {
        input:
        in_vcf_list=VCF_LIST,
        in_cores=CORES,
        in_disk=DISK,
        in_mem=MEM,
        in_preemptible=PREEMPTIBLE
    }
    
    call normalize {
        input:
        in_vcf_file=merge.out_vcf,
        in_fasta_file=REF_FASTA,
        in_fai_file=REF_FASTA_IDX,
        in_cores=CORES,
        in_disk=DISK,
        in_mem=MEM,
        in_preemptible=PREEMPTIBLE
    }

    call sort {
        input:
        in_vcf_file=normalize.out_vcf,
        in_cores=CORES,
        in_disk=DISK,
        in_mem=MEM,
        in_preemptible=PREEMPTIBLE
    }

    output {
        File output_vcf = sort.out_vcf
        File output_vcf_idx = sort.out_vcf_idx
    }
}

task merge {
    input {
        Array[File] in_vcf_list
        Int in_cores
        Int in_disk
        Int in_mem
        Int in_preemptible
    }
    command <<<
    set -eux -o pipefail

    rm -f vcf_list.txt
    while read invcf
    do
        outvcf=`basename $invcf`
        bcftools norm -m -both $invcf | bcftools view --exclude 'GT="0/0" || GT="0" || GT~"\."' -O z > $outvcf
        bcftools index $outvcf
        echo $outvcf >> vcf_list.txt
    done < ~{write_lines(in_vcf_list)}

    bcftools merge -0 --threads ~{in_cores} -m none -l vcf_list.txt -O z > merged.vcf.bgz
    >>>
    output {
        File out_vcf= "merged.vcf.bgz"
    }
    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10--h5d15f04_0"
        memory: in_mem + " GB"
        cpu: 1
        disks: "local-disk " + in_disk + " SSD"
        preemptible: in_preemptible
    }
}

task normalize {
    input {
        File in_vcf_file
        File in_fasta_file
        File in_fai_file
        Int in_cores
        Int in_disk
        Int in_mem
        Int in_preemptible
    }
    command <<<
    set -eux -o pipefail
    bcftools norm --threads ~{in_cores} -f ~{in_fasta_file} -O z ~{in_vcf_file} > norm.vcf.bgz
    >>>
    output {
        File out_vcf= "norm.vcf.bgz"
    }
    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10--h5d15f04_0"
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
        preemptible: in_preemptible
    }
}

task sort {
    input {
        File in_vcf_file
        Int in_cores
        Int in_disk
        Int in_mem
        Int in_preemptible
    }
    command <<<
    set -eux -o pipefail

    mkdir -p temp
    bcftools sort -T temp -m 2G -O z  ~{in_vcf_file} > sorted.vcf.bgz
    bcftools index -t sorted.vcf.bgz
    >>>
    output {
        File out_vcf= "sorted.vcf.bgz"
        File out_vcf_idx= "sorted.vcf.bgz.tbi"
    }
    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10--h5d15f04_0"
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
        preemptible: in_preemptible
    }
}

