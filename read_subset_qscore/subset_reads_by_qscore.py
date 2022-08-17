import argparse
import gzip


parser = argparse.ArgumentParser(description='Subset reads by qscore.')
parser.add_argument('-s', help='gzipped summary file from Guppy',
                    required=True)
parser.add_argument('-f', help='gzipped fastq file from Guppy', required=True)
parser.add_argument('-q', default=10, type=int, help='minimum qscore')

args = parser.parse_args()

# read summary file and save reads that should be filtered out
# (hopefully less than reads to keep)
sumf = gzip.open(args.s, 'rt')
reads_to_rm = set()
# read headers to know which columns have read name and qscore
sum_head = sumf.readline().rstrip().split('\t')
q_col_idx = sum_head.index('mean_qscore_template')
rn_col_idx = sum_head.index('read_id')
for line in sumf:
    line = line.rstrip().split('\t')
    if float(line[q_col_idx]) < args.q:
        reads_to_rm.add(line[rn_col_idx])
sumf.close()

# stream through the input fastq
infq = gzip.open(args.f, 'rt')
fq_l = 3         # which of the 4 FASTQ lines are we at
skip = False     # are we skipping this read
for infq_line in infq:
    # if not the read name line, skip or write info
    if fq_l < 3:
        if not skip:
            print(infq_line.rstrip())
        fq_l += 1
        continue
    # else it's the read name line
    fq_l = 0
    # get read name
    rname = infq_line.split()[0][1:]
    # skip or not
    if rname in reads_to_rm:
        skip = True
    else:
        skip = False
        print(infq_line.rstrip())

# close file connection
infq.close()
