import sys
import os
import gzip
import argparse
import pandas as pd

infile = ""
outpath = ""
parser = argparse.ArgumentParser()
parser.add_argument('-vcf','--vcf', help='Path to the gnomAD_SV vcf file', required=True)
parser.add_argument('-o','--out', help='Path to the output csv file', required=True)
args = vars(parser.parse_args())
if "vcf.gz.vcf" in args['vcf']:
	infile = args['vcf']
	a = os.path.splitext(infile)
	infile = a[0]
if ".csv" in args['out']:
	outpath = args['out']
else:
	outpath = args['out']+".csv"
print(infile)
counter = 0
flag = False
chroms = []
pos = []
variant_ids = []
AFs = []
with gzip.open(infile, 'rb') as f:
    for line in f:
        # print(line.decode().strip())
        l = line.decode().strip()
        if "#CHROM" in line.decode().strip():
        	flag = True
        	continue
        if flag:
        	a = l.strip().split(';')
        	aa = a[0].split('\t')
        	chroms.append(str(aa[0]))
        	pos.append(str(aa[1]))
        	variant_ids.append(str(aa[2]))
        	for i in a:
        		if "AF=" in i and not "_AF=" in i:
        			AFs.append(str(i.split("=")[-1]))
dict_ = {'Chrom':chroms, "Pos":pos, "ID":variant_ids, "AF":AFs}
df = pd.DataFrame(dict_)
df = df.reindex(columns=['Chrom', 'Pos', 'ID', 'AF'])
df.to_csv(outpath, index=False)


