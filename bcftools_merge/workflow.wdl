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
    }

    call merge {
        input:
        in_vcf_list=VCF_LIST,
        in_disk=DISK,
        in_mem=MEM
    }
    
    call normalize {
        input:
        in_vcf_file=merge.out_vcf,
        in_fasta_file=REF_FASTA,
        in_fai_file=REF_FASTA_IDX,
        in_cores=CORES,
        in_disk=DISK,
        in_mem=MEM
    }
    
    output {
        File output_vcf = normalize.out_vcf
        File output_vcf_idx = normalize.out_vcf_idx
    }
}

task merge {
    input {
        Array[File] in_vcf_list
        Int in_disk
        Int in_mem
    }
    command <<<
    set -eux -o pipefail
    
    rm -f vcf.list.txt
    for INVCF in ~{sep=" " in_vcf_list}
    do
        echo "Indexing $INVCF..."
        bcftools index $INVCF
        echo $INVCF >> vcf.list.txt
    done

    echo "Merge VCFs and sort..."
    bcftools merge -m all -l vcf.list.txt | bcftools sort -O z > merged.vcf.bgz
    >>>
    output {
        File out_vcf= "merged.vcf.bgz"
    }
    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10--h5d15f04_0"
        memory: in_mem + " GB"
        cpu: 1
        disks: "local-disk " + in_disk + " SSD"
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
    }
    command <<<
    set -eux -o pipefail

    bcftools norm -m -both --threads ~{in_cores} -f ~{in_fasta_file} ~{in_vcf_file} | bcftools sort -O z > norm.vcf.bgz
    bcftools index -t norm.vcf.bgz
    >>>
    output {
        File out_vcf= "norm.vcf.bgz"
        File out_vcf_idx= "norm.vcf.bgz.tbi"
    }
    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10--h5d15f04_0"
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
    }
}

