Subset a fastq for reads with a minimum basecalling qscore. 

Both fastq and summary inputs should have the reads in the same order. They should if were the output of the same Guppy run. 

## Testing locally

Test with :

- [test.reads.fastq.gz](../read_stats/test.reads.fastq.gz) from the read_stats workflow's test data
- [test.summary.tsv.gz](test.summary.tsv.gz)

~~~
java -jar $CROMWELL_JAR run workflow.wdl -i inputs.json
~~~

## Docker container

The workflow uses the [quay.io/jmonlong/read_subset_qscore](https://quay.io/repository/jmonlong/read_subset_qscore).
It contains python3 and the [subset_reads_by_qscore.py](subset_reads_by_qscore.py) script.
It's made with the [Dockerfile](Dockerfile) in this repo.

