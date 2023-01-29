## Test locally

```
## download a small VCF to test (slice of HG002 GIAB)
bcftools view ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz chr7:97352092-104568891 > test.vcf
bgzip test.vcf

## download GRCh38.105 database
wget https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_GRCh38.105.zip

java -jar $CROMWELL_JAR run workflow.wdl -i test.inputs.json
```
