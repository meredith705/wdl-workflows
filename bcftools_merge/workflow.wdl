version 1.0

workflow bcftools_merge {
    meta {
	author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Merge VCFs using bcftools. Multi-allelic records are split and normalized (requires a reference fasta) before being merged."
    }
    input {
        Array[File] VCF_LIST
        File REF_FASTA
        File REF_FASTA_IDX
        Boolean SPLIT_MODE = true  ## split/normalized/realigned variants? If not, only (bcftools) merge
        Int CORES = 1
        Int DISK = 100
        Int MEM = 10
        Int PREEMPTIBLE = 0
    }

    call merge {
        input:
        in_vcf_list=VCF_LIST,
        in_fasta_file=REF_FASTA,
        in_fai_file=REF_FASTA_IDX,
        in_split_mode=SPLIT_MODE,
        in_cores=CORES,
        in_disk=DISK,
        in_mem=MEM,
        in_preemptible=PREEMPTIBLE
    }
    
    output {
        File output_vcf = merge.out_vcf
        File output_vcf_idx = merge.out_vcf_idx
    }
}

task merge {
    input {
        Array[File] in_vcf_list
        File in_fasta_file
        File in_fai_file
        Boolean in_split_mode
        Int in_cores
        Int in_disk
        Int in_mem
        Int in_preemptible
    }
    command <<<
    set -eux -o pipefail

    # prepare scripts to split and normalize SVs for each sample
    rm -f vcf_list.txt
    mkdir -p temp scripts
    while read invcf
    do
        outvcf=`basename $invcf`
        if [ ~{in_split_mode} == true ]; then
            # split variants
            echo "bcftools norm -m -both $invcf -O z > $outvcf.split.vcf.gz" > scripts/script_samp.$outvcf.sh
            # realign and split variants
            echo "python3 /scripts/align_variants.py -i $outvcf.split.vcf.gz -f ~{in_fasta_file} -o $outvcf.al.vcf" >> scripts/script_samp.$outvcf.sh
            # filter ref calls, normalize, sort and merge hets that match exactly into homs
            echo "bcftools view --exclude 'GT=\"0/0\" || GT=\"0\" || GT~\"\\.\"' $outvcf.al.vcf | bcftools norm -f ~{in_fasta_file} | bcftools sort -T temp/$outvcf.temp | python3 /scripts/merge_exact_hets.py | bgzip > $outvcf" >> scripts/script_samp.$outvcf.sh
            # index new VCF
            echo "tabix -f $outvcf" >> scripts/script_samp.$outvcf.sh
            # clean up
            echo "rm -f $outvcf.split.vcf.gz $outvcf.al.vcf" >> scripts/script_samp.$outvcf.sh
            echo $outvcf >> vcf_list.txt
        else
            # no splitting/normalization/realignment, just index input VCF
            echo "bcftools index $invcf" > scripts/script_samp.$outvcf.sh
            echo $invcf >> vcf_list.txt
        fi
    done < ~{write_lines(in_vcf_list)}

    # run the scripts in parallel
    ls scripts/script_samp.*.sh | parallel -j ~{in_cores} sh {} 

    # merge the SVs across samples
    bcftools merge -0 --threads ~{in_cores} -m none -l vcf_list.txt -O z > merged.vcf.bgz
    tabix merged.vcf.bgz
    >>>
    output {
        File out_vcf= "merged.vcf.bgz"
        File out_vcf_idx= "merged.vcf.bgz.tbi"
    }
    runtime {
        docker: "jmonlong/merge-sv-vg@sha256:f10085ba57c7f5f26e6bb066e5cd3befe319ece6d547d64c81224abfc8a8b1d4"
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
        preemptible: in_preemptible
    }
}
