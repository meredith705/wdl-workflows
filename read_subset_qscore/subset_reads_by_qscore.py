import argparse
import gzip
import datetime
import tracemalloc


# function to print the progress
def printProgress(reads_s, reads_t):
    print('{} - {} reads skipped out of {} ({}%)'.format(
        datetime.datetime.now(),
        reads_s,
        reads_t,
        round(100 * reads_s / reads_t, 2)))


parser = argparse.ArgumentParser(description='Subset reads by qscore.')
parser.add_argument('-s', help='gzipped summary file from Guppy',
                    required=True)
parser.add_argument('-f', help='gzipped fastq file from Guppy', required=True)
parser.add_argument('-q', default=10, type=int, help='minimum qscore')
parser.add_argument('-o', help='the output gzipped fastq file')

args = parser.parse_args()
infq = gzip.open(args.f, 'rt')
outfq = gzip.open(args.o, 'wt')

# read summary file and save reads that should be filtered out
# (hopefully less than reads to keep)
sumf = gzip.open(args.s, 'rt')
reads_to_rm = set()
# measure the memory footprint of saving those read names
tracemalloc.start()
## read headers to know which columns have read name and qscore
sum_head = sumf.readline().rstrip().split('\t')
q_col_idx = sum_head.index('mean_qscore_template')
rn_col_idx = sum_head.index('read_id')
for line in sumf:
    line = line.rstrip().split('\t')
    if float(line[q_col_idx]) < args.q:
        reads_to_rm.add(line[rn_col_idx])
sumf.close()

print('{} - {} reads to remove (if in the fastq)'.format(
    datetime.datetime.now(),
    len(reads_to_rm)))

# print the peak memory usage so far
mem_size, peak_mem_size = tracemalloc.get_traced_memory()
mem_gb = round(peak_mem_size/(1024*1024*1024), 5)
print('{}Gb of memory used to store the read names'.format(mem_gb))
# stopping memory monitoring
tracemalloc.stop()

# save some stats to report at the end
total_reads = 0
total_reads_skipped = 0

# stream through the input fastq
fq_l = 3         # which of the 4 FASTQ lines are we at
skip = False     # are we skipping this read
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
    # skip or not
    if rname in reads_to_rm:
        skip = True
        total_reads_skipped += 1
    else:
        skip = False
        outfq.write(infq_line)
    total_reads += 1
    if total_reads % 500000 == 0:
        printProgress(total_reads_skipped, total_reads)

# close file connections
infq.close()
outfq.close()

# final stats
printProgress(total_reads_skipped, total_reads)
