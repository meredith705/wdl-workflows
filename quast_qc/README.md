Run QUAST to get QC metrics on an assembly slightly tweaked from [human-pangenomics/hpp_production_workflows/Quast](https://dockstore.org/workflows/github.com/human-pangenomics/hpp_production_workflows/Quast:master?tab=info) to allow for an input array of fasta files (that will be concatenated).

Test locally with

~~~
java -jar $CROMWELL_JAR run workflow.wdl -i inputs.json
~~~
