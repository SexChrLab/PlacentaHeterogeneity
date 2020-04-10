# This script calculates the allele balance
import argparse

parser = argparse.ArgumentParser(description='Calculate allele balance (unphased)')
parser.add_argument('--input',required=True,help='Input the path to the output from asereadcounter')
parser.add_argument('--output',required=True,help='Input the path to the output file')

args = parser.parse_args()

# Header of input file: contig  position        variantID       refAllele       altAllele       refCount        altCount        totalCount      lowMAPQDepth    lowBaseQDepth   rawDepth        otherBases      improperPairs
input = args.input

# Format of output file:
output = open(args.output, 'w')
header = ['chr', 'position', 'ref_allele', 'alt_allele', 'ref_count', 'alt_count', 'total_count', 'allele_balance']
print ('\t'.join(header), file=output)

with open(input, 'r') as f:
    for line in f:
        if not line.startswith('contig'):
            items = line.rstrip('\n').split('\t')
            out = [items[0], items[1], items[3], items[4], items[5], items[6], items[7]]
            ref_ratio = float(items[5])/float(items[7])
            alt_ratio = float(items[6])/float(items[7])

            if ref_ratio > alt_ratio:
                out.append(str(ref_ratio))
            else:
                out.append(str(alt_ratio))
            print ('\t'.join(out), file=output)

