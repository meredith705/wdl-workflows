version 1.0

workflow minimap2 {
    meta {
	author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Align two sets of sequences using minimap2 and return a gzipped .paf file. For example, an assembly vs a reference genome." 
    }
    input {
        File TARGET_SEQ_FILE
        File QUERY_SEQ_FILE
        String MINIMAP2_CONTAINER = "quay.io/biocontainers/minimap2:2.17--hed695b0_3"
        String MINIMAP2_ARGS=""
        String OUT_LABEL="minimap2-out"
        Int CORES = 16
        Int DISK = 500
        Int MEM = 50
    }

    call alignAndGzip {
        input:
        in_target_seq_file=TARGET_SEQ_FILE,
        in_query_seq_file=QUERY_SEQ_FILE,
        in_container=MINIMAP2_CONTAINER,
        in_args=MINIMAP2_ARGS,
        in_out_label=OUT_LABEL,
        in_cores=CORES,
        in_disk=DISK,
        in_mem=MEM
    }
    
    output {
        File output_paf = alignAndGzip.out_paf
    }
}

task alignAndGzip {
    input {
        File in_target_seq_file
        File in_query_seq_file
        String in_container
        String in_args
        String in_out_label
        Int in_cores
        Int in_disk
        Int in_mem
    }
    command <<<
    set -eux -o pipefail
    minimap2 -t ~{in_cores} ~{in_args} ~{in_target_seq_file} ~{in_query_seq_file} | gzip > ~{in_out_label}.paf.gz
    >>>
    output {
        File out_paf= "~{in_out_label}.paf.gz"
    }
    runtime {
        docker: in_container
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
    }
}

