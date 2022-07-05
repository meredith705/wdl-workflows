Uses [Liger2LiGer](https://github.com/rlorigro/Liger2LiGer) to estimate the amount of chimeric reads in a set of long reads in FASTQ format.

- Docker image deposited at [quay.io/jmonlong/liger2liger](https://quay.io/repository/jmonlong/liger2liger) ([Dockerfile](https://github.com/jmonlong/docker-liger2liger))
- Test data:
  - [ref.fa.gz](ref.fa.gz) fasta with 2 reference contigs
  - [test.reads.fastq.gz](test.reads.fastq.gz) FASTQ with two reads, one chimeric, one non-chimeric.

~~~
java -jar $CROMWELL_JAR run workflow.wdl -i inputs.json
~~~
