library(dplyr)

nsamps = 1000
nvars = 1000
outf = 'test_input.chr1.vcf'
pos.max = 248956422

## headers
cat('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n', file=outf)
hh = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', paste0('s', 1:nsamps))
cat(paste0(paste(hh, collapse='\t'), '\n'), file=outf, append=TRUE)

## variants
pos = sort(unique(round(runif(nvars, 0, pos.max))))
rec = tibble(chr='chr1', pos=pos, id='.', r='A', a='T', q='.', f='.', i='.', format='GT')
gts = matrix(sample(c('0/0', '0/1', '1/1'), length(pos)*nsamps, replace=TRUE, prob=c(.5,.4,.1)), length(pos))
rec = cbind(rec, gts)
write.table(rec, file=outf, append=TRUE, sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
