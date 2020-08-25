wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
zcat ce11.fa.gz | head -10000 | gzip > target.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/ce2/bigZips/ce2.fa.gz
zcat ce2.fa.gz | head -10000 | gzip > query.fa.gz

java -jar ~/soft/cromwell-43.jar run workflow.wdl -i inputs.json
