import argparse
import gzip

parser = argparse.ArgumentParser(description='Subset reads by qscore. ' +
                                 'Reads should be in the same order in' +
                                 ' the summary and input fastq files.')
parser.add_argument('-s', help='gzipped summary file from Guppy',
                    required=True)
parser.add_argument('-f', help='gzipped fastq file from Guppy', required=True)
parser.add_argument('-q', default=10, type=int, help='minimum qscore')
parser.add_argument('-o', help='the output gzipped fastq file')

args = parser.parse_args()
sumf = gzip.open(args.s, 'rt')
infq = gzip.open(args.f, 'rt')
outfq = gzip.open(args.o, 'wt')

# read header of summary file and find column id with qscore
sum_head = sumf.readline().rstrip().split('\t')
q_col_idx = sum_head.index('mean_qscore_template')
rn_col_idx = sum_head.index('read_id')

# save some stats to report at the end
total_reads = 0
total_reads_skipped = 0

# stream through the input fastq
fq_l = 3         # which of the 4 FASTQ lines are we at
skip = False     # are we skipping this read
cur_sum_rn = ''  # current read name in the summary file
for infq_line in infq:
    # if not the read name line, skip or write info
    if fq_l < 3:
        if not skip:
            outfq.write(infq_line)
        fq_l += 1
        continue
    # else it's the read name line
    fq_l = 0
    # get read name
    rname = infq_line.split()[0][1:]
    # look for the read in the summary file
    while cur_sum_rn != rname:
        sum_line = sumf.readline().rstrip().split('\t')
        cur_sum_rn = sum_line[rn_col_idx]
    # check qscore vs min threshold
    if float(sum_line[q_col_idx]) >= args.q:
        skip = False
        outfq.write(infq_line)
    else:
        skip = True
        total_reads_skipped += 1
    total_reads += 1

# close file connections
sumf.close()
infq.close()
outfq.close()

# print some stats
print('{} reads skipped out of {} ({}%)'.format(
    total_reads_skipped,
    total_reads,
    round(100 * total_reads_skipped / total_reads, 2)))
