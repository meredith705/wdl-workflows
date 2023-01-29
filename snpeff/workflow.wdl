version 1.0

workflow snpeff_annotate {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Annotate a VCF with SNPeff"
    }
    input {
        File VCF
        File SNPEFF_DB
        File SNPEFF_DB_NAME
        Boolean SPLIT_MULTIAL = false
        Boolean SORT_INDEX_VCF = false
        Boolean KEEP_ANN_ONLY = true
    }

    if (SPLIT_MULTIAL){
        call split_multiallelic_vcf {
            input:
            input_vcf=VCF
        }
    }

    File current_vcf = select_first([split_multiallelic_vcf.vcf, VCF])
    
    call annotate_with_snpeff {
        input:
        input_vcf=current_vcf,
        snpeff_db=SNPEFF_DB,
        db_name=SNPEFF_DB_NAME
    }

    if (KEEP_ANN_ONLY || SORT_INDEX_VCF){
        call format_vcf {
            input:
            input_vcf=annotate_with_snpeff.vcf,
            only_ann=KEEP_ANN_ONLY,
            sort_index=SORT_INDEX_VCF
        }
    }
    
    File final_vcf = select_first([format_vcf.vcf, annotate_with_snpeff.vcf])
    
    output {
        File vcf = final_vcf
        File? vcf_index = format_vcf.vcf_index
    }
}

task split_multiallelic_vcf {
    input {
        File input_vcf
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(input_vcf, "GB")) + 20
    }

    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
    command <<<
    set -eux -o pipefail
    
    bcftools norm -m -both --threads ~{threadCount} -Oz -o ~{basen}.norm.vcf.gz ~{input_vcf}
    >>>
    
    output {
        File vcf = "~{basen}.norm.vcf.gz"
    }
    
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.16--hfe4b78e_1"
        preemptible: 1
    }
}

task format_vcf {
    input {
        File input_vcf
        Boolean only_ann
        Boolean sort_index
        Int memSizeGB = 4
        Int diskSizeGB = 5*round(size(input_vcf, "GB")) + 20
    }

    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
    command <<<
    set -eux -o pipefail

    VCF=~{input_vcf}
    if [ ~{only_ann} == true ]
    then
        bcftools view -i "COUNT(INFO/ANN)>0" -Oz -o ~{basen}.formatted.vcf.gz ~{input_vcf}
        VCF=~{basen}.formatted.vcf.gz
    fi

    if [ ~{only_ann} == true ]
    then
        cp $VCF temp_$VCF
        bcftools sort -Oz -o ~{basen}.formatted.vcf.gz temp_$VCF
        bcftools index -t -o ~{basen}.formatted.vcf.gz.tbi ~{basen}.formatted.vcf.gz
    fi
    >>>
    
    output {
        File vcf = "~{basen}.formatted.vcf.gz"
        File? vcf_index = "~{basen}.formatted.vcf.gz.tbi"
    }
    
    runtime {
        memory: memSizeGB + " GB"
        cpu: 1
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/bcftools:1.16--hfe4b78e_1"
        preemptible: 1
    }
}

task annotate_with_snpeff {
    input {
        File input_vcf
        File snpeff_db
        String db_name
        Int memSizeGB = 16
        Int threadCount = 2
        Int diskSizeGB = 5*round(size(input_vcf, "GB") + size(snpeff_db, 'GB')) + 20
    }

    Int snpeffMem = if memSizeGB < 6 then 2 else memSizeGB - 4
    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
	command <<<
        set -eux -o pipefail

        unzip ~{snpeff_db}
        
        snpEff -Xmx~{snpeffMem}g -nodownload -no-intergenic \
               -dataDir ${PWD}/data ~{db_name} \
               ~{input_vcf} | gzip > ~{basen}.snpeff.vcf.gz
	>>>

	output {
		File vcf = "~{basen}.snpeff.vcf.gz"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/snpeff:5.1d--hdfd78af_0"
        preemptible: 1
    }
}
