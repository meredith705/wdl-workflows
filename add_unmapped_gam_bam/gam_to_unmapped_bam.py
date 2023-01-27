import argparse
import fileinput


parser = argparse.ArgumentParser(description='Extract reads from GAM/json' +
                                 'and output as unmapped SAM.')
parser.add_argument('-r', help='the read names', required=True)
args = parser.parse_args()

# read names
rnames = {1}
with open(args.r, 'rt') as inf:
    for line in inf:
        rnames.add(line.rstrip())

for line in fileinput.input(files='-'):
    line = line.rstrip().split('\t')
    if line[0] in rnames:
        print('{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}'.format(line[0], line[1],
                                                          "~"*len(line[1])))
