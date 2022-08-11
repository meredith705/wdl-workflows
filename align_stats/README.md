Uses [wambam](https://github.com/rlorigro/wambam) to compute the distribution of alignment identity for reads.
This workflow aligns reads with [minimap2](https://github.com/lh3/minimap2) and then run wambam.

- Docker image deposited at [quay.io/jmonlong/wambam](https://quay.io/repository/jmonlong/wambam) ([Dockerfile](https://github.com/jmonlong/docker-wambam))
- Test data from the [*chimera_stats*](../chimera_stats) workflow:
  - [ref.fa.gz](../chimera_stats/ref.fa.gz) fasta with 2 reference contigs
  - [test.reads.fastq.gz](../chimera_stats/test.reads.fastq.gz) FASTQ with two reads, one chimeric, one non-chimeric.

~~~
java -jar $CROMWELL_JAR run workflow.wdl -i inputs.json
~~~
