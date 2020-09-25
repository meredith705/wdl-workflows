ln -s ~/Documents/datagenome/hg38.fa
ln -s ~/Documents/datagenome/hg38.fa.fai 

java -jar ~/soft/cromwell-43.jar run workflow.wdl -i inputs.json
