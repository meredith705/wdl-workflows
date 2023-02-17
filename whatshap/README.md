## Phasing evaluation with WhatsHap

Simple workflow to evaluate the phasing in a VCF based on a phased truthset VCF, using WhatsHap compare.

## Docker container

```
docker build -t whatshap .
docker tag whatshap quay.io/jmonlong/whatshap:1.7
docker push quay.io/jmonlong/whatshap:1.7
```

## Test locally

```
## download a small VCF to test (slice of HG002 GIAB)
bcftools view https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.vcf.gz chr7:97352092-104568891 > test.vcf
bgzip test.vcf

miniwdl run --as-me --copy-input-files workflow.wdl -i test.inputs.json
```
