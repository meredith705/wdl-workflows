version 1.0

workflow WhatsHapCompare {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Evaluate phasing in VCF using WhatsHap. Adapted from Trevor Pesout's version."
    }
    input {
        File CALL_VCF
        File TRUTH_VCF
        String? SAMPLE_NAME = "sample"
    }

    call whatshap_compare {
        input:
        call_vcf=CALL_VCF,
        truth_vcf=TRUTH_VCF,
        sample_name=SAMPLE_NAME
    }

    output {
        File out_tarball = whatshap_compare.tarball
        File out_summary = whatshap_compare.summary
        File out_stats = whatshap_compare.stats
    }
}

task whatshap_compare {
    input {
        File call_vcf
        File truth_vcf
        String? sample_name = "sample"
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(call_vcf, "GB") + size(truth_vcf, "GB")) + 20
    }

    String call_pref = sub(sub(basename(call_vcf), ".gz", ""), ".vcf", "")
    String truth_pref = sub(sub(basename(truth_vcf), ".gz", ""), ".vcf", "")
    command <<<
    set -eux -o pipefail

    echo ~{sample_name} > sample.txt
    
    ## format sample name for truth VCF
    bcftools reheader -s sample.txt ~{truth_vcf} > ~{truth_pref}.vcf
   
    ## format sample name for call VCF
    bcftools reheader -s sample.txt ~{call_vcf} > ~{call_pref}.vcf

    mkdir whatshap_eval_results
    
    # get stats
    whatshap stats \
             --tsv whatshap_eval_results/~{sample_name}.stats.tsv \
             --block-list whatshap_eval_results/~{sample_name}.blocks.txt \
             --chr-lengths /opt/whatshap/chr_lengths \
             ~{call_pref}.vcf \
             > whatshap_eval_results/~{sample_name}.stats.txt
    
    # get comparison
    whatshap compare \
             --tsv-pairwise whatshap_eval_results/~{sample_name}.pairwise.tsv \
             --tsv-multiway whatshap_eval_results/~{sample_name}.multiway.tsv \
             --switch-error-bed whatshap_eval_results/~{sample_name}.switch_error.bed \
             --longest-block-tsv whatshap_eval_results/~{sample_name}.longest_block.tsv \
             ~{call_pref}.vcf ~{truth_pref}.vcf \
             > whatshap_eval_results/~{sample_name}.compare.txt
    
    # tarball it
    tar czvf whatshap_eval_results.~{sample_name}.tar.gz whatshap_eval_results/

    >>>
    
    output {
        File tarball = "whatshap_eval_results.~{sample_name}.tar.gz"
        File stats = "whatshap_eval_results/~{sample_name}.stats.txt"
        File summary = "whatshap_eval_results/~{sample_name}.pairwise.tsv"
    }
    
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/whatshap@sha256:03e5d1126b46c9bc671219b4eb7e9647fcb8daf22ce68acaac7496a49d8e057d"
        preemptible: 1
    }
}
