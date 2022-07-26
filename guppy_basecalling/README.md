Basecall a fast5 file using Guppy. 

This workflow was adapted from the [scatterGuppy.wdl in the human-pangenomics/hpp_production_workflows repo](https://github.com/human-pangenomics/hpp_production_workflows/blob/master/data_processing/wdl/scatterGuppy.wdl).
It was written originally by Jimin Park.

I think Jimin created the docker container from [this Dockerfile](https://github.com/human-pangenomics/hpp_production_workflows/blob/master/data_processing/docker/guppy/Dockerfile).

My tweaks to the workflow:
- Gzip the output "read summary" (takes about 2 Gb otherwise)
- Make it a simpler single-sample workflow, i.e. one input fast5.tar, one output bam, etc

The small fast5 file used to test the workflow was downloaded from [https://github.com/mbhall88/fast5seek/tree/master/tests/data/fast5](https://github.com/mbhall88/fast5seek/tree/master/tests/data/fast5).

~~~
java -jar $CROMWELL_JAR run workflow.wdl -i inputs.json
~~~
