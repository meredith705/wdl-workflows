import argparse
from cyvcf2 import VCF, Writer


parser = argparse.ArgumentParser()
parser.add_argument('-v', help='variants in VCF (can be bgzipped)',
                    required=True)
parser.add_argument('-o', default='vcf_stats.tsv',
                    help='output TSV with stats')
args = parser.parse_args()

del_het = del_hom = ins_het = ins_hom = 0
ps_coord = {}

# Read VCF
vcf = VCF(args.v)
for variant in vcf:
    # filter homozygous ref variants
    if variant.num_het == 0 and variant.num_hom_alt == 0:
        continue
    # update the boundary of phase block is PS is present
    if 'PS' in variant.FORMAT:
        ps = variant.format('PS')[0][0]
        if ps in ps_coord:
            ps_coord[ps]['start'] = min(variant.start, ps_coord[ps]['start'])
            ps_coord[ps]['end'] = max(variant.end, ps_coord[ps]['end'])
        else:
            ps_coord[ps] = {'start': variant.start, 'end': variant.end}
    # if SV, count if heterozygous or homozygous
    if len(variant.REF) > 30 or len(variant.ALT[0]) > 30:
        if len(variant.REF) > len(variant.ALT[0]):
            # deletion
            del_het += variant.num_het
            del_hom += variant.num_hom_alt
        else:
            # insertion    
            ins_het += variant.num_het
            ins_hom += variant.num_hom_alt
vcf.close()

# compute phase block stats
ps_sizes = []
for ps in ps_coord:
    ps_sizes.append(ps_coord[ps]['end'] - ps_coord[ps]['start'])
ps_sizes = sorted(ps_sizes)
half_tot_size = sum(ps_sizes)/2
ps_n50 = 0
ps_csum = 0
for pss in ps_sizes:
    if ps_csum > half_tot_size:
        break
    ps_n50 = pss
    ps_csum += pss

# write output TSV
with open(args.o, 'wt') as outf:
    headers = ['phase_block_n50', 'del_het', 'del_hom', 'ins_het', 'ins_hom']
    outf.write("\t".join(headers) + '\n')
    outf.write("{}\t{}\t{}\t{}\t{}\n".format(ps_n50, del_het, del_hom, ins_het, ins_hom))
