# In this script, we are going to subset the file /scratch/tphung3/PlacentaSexDiff/A_placenta/03_analyze_ase/results/chrX/OBG0178_OBG0178-1-020_chrX_allele_balance.tsv to contain only positions in xist
import sys
import os
xist_start = 73820892
xist_end = 73851867

for fn in os.listdir('/scratch/tphung3/PlacentaSexDiff/A_placenta/02_run_asereadcounter/asereadcounter/HISAT/chrX/'):
    if fn.endswith('_chrX.tsv'):
        print (fn)
        with open(os.path.join('/scratch/tphung3/PlacentaSexDiff/A_placenta/02_run_asereadcounter/asereadcounter/HISAT/chrX/', fn), 'r') as f:
            for line in f:
                if line.startswith('chrX'):
                    items = line.rstrip('\n').split('\t')
                    if int(items[1]) >= xist_start and int(items[1]) <= xist_end:
                        print (fn)
                        print (line.rstrip('\n'))