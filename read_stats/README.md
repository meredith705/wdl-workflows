Simple stats about reads in one gzipped FASTQ file. Simplified from the [*read_stats* HPP workflow](https://github.com/human-pangenomics/hpp_production_workflows/blob/master/QC/wdl/tasks/read_stats.wdl) (one file input, must be gzipped FASTQ, but all in one task).

Test with [test.reads.fastq.gz](test.reads.fastq.gz):

~~~
java -jar $CROMWELL_JAR run workflow.wdl -i inputs.json
~~~
