version 1.0

workflow beagle_phasing {
    meta {
	author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Phasing a VCF using BEAGLE 5.1"
    }
    input {
        File IN_VCF
        String IN_REF_VERSION = "GRCh38" # options: GRCh36, GRCh37, GRCh38
        String IN_CHR
        Int CORES = 1
        Int DISK = 100
        Int MEM = 10
        Int PREEMPTIBLE = 0
    }

    call phasing {
        input:
        in_vcf=IN_VCF,
        in_ref_version=IN_REF_VERSION,
        in_chr=IN_CHR,
        in_cores=CORES,
        in_disk=DISK,
        in_mem=MEM,
        in_preemptible=PREEMPTIBLE
    }
    
    output {
        File output_vcf = phasing.out_vcf
    }
}

task phasing {
    input {
        File in_vcf
        String in_ref_version
        String in_chr
        Int in_cores
        Int in_disk
        Int in_mem
        Int in_preemptible
    }
    command <<<
    set -eux -o pipefail

    ## download genetic map for specified reference version
    wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.~{in_ref_version}.map.zip
    unzip plink.~{in_ref_version}.map.zip

    ## add 'chr' prefix to genetic map
    awk '{print "chr"$0}' plink.~{in_chr}.~{in_ref_version}.map > genetic.map

    ## run beagle
    java -Xmx~{in_mem}g -jar /beagle/beagle.jar gt=~{in_vcf} out=phased map=genetic.map nthreads=~{in_cores}
    >>>
    output {
        File out_vcf= "phased.vcf.gz"
    }
    runtime {
        docker: "jmonlong/beagle"
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
        preemptible: in_preemptible
    }
}
