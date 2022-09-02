Uses [Liger2LiGer](https://github.com/rlorigro/Liger2LiGer) to estimate the amount of chimeric reads in a set of long reads in FASTQ format or BAM.

- Docker image deposited at [quay.io/jmonlong/liger2liger](https://quay.io/repository/jmonlong/liger2liger) ([Dockerfile](https://github.com/jmonlong/docker-liger2liger))
- Inputs are either:
  - a `FASTQ_FILE` and `REFERENCE_FILE` reference fasta file (reads will be aligned to the reference with minimap2)
  - a `BAM_FILE` BAM file with the aligned reads.
- Test data:
  - when input is reads + reference fasta
      - [ref.fa.gz](ref.fa.gz) fasta with 2 reference contigs
      - [test.reads.fastq.gz](test.reads.fastq.gz) FASTQ with two reads, one chimeric, one non-chimeric.
  - when input is a BAM file
      - [test.reads.bam](test.reads.bam) BAM with the same two reads, one chimeric, one non-chimeric.
      
~~~
java -jar $CROMWELL_JAR run workflow.wdl -i inputs.json
## or
java -jar $CROMWELL_JAR run workflow.wdl -i inputs.bam.json
~~~
