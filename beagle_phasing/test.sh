Rscript simvcf.R
bgzip -f test_input.chr1.vcf

wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip plink.GRCh38.map.zip

java -Xmx4g -jar /beagle/beagle.jar gt=test_input.chr1.vcf.gz out=phased map=genetic.map nthreads=1
