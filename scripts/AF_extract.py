import sys
import os
import gzip
import argparse
import pandas as pd
import re

infile = ""
infile_ = ""
outpath = ""
parser = argparse.ArgumentParser()
parser.add_argument('-vcf','--vcf', help='Path to the gnomAD_SV vcf file', required=True)
parser.add_argument('-out','--out', help='Path to the output tsv file', required=True)
args = vars(parser.parse_args())
if "vcf.gz.vcf" in args['vcf']:
	infile = args['vcf']
	a = os.path.splitext(infile)
	infile_ = a[0]
	os.system("cp "+infile+" "+infile_)
if (args['out']).endswith(".tsv"):
	outpath = args['out']
else:
	outpath = args['out']+".tsv"
print(infile_)
counter = 0
flag = False
chroms = []
pos = []
variant_ids = []
AFs = []
types = []
ends = []
with gzip.open(infile_, 'rb') as f:
    for line in f:
        # print(line.decode().strip())
        l = line.decode().strip()
        if "#CHROM" in line.decode().strip():
        	flag = True
        	continue
        if flag:
        	a = re.split(';|\t', l)
        	# aa = a[0].strip().split('\t')
        	chroms.append("chr"+str(a[0]))
        	pos.append(str(a[1]))
        	variant_ids.append(str(a[2]))
        	types.append(str(a[4]))
        	for i in a:
        		if "END=" in i:
        			ends.append(str(i.split('=')[-1]))
        	# ends.append(str(a[7].split('=')[-1]))
        	for i in a:
        		if "AF=" in i and not "_AF=" in i:
        			AFs.append(str(i.split("=")[-1]))
dict_ = {'chr':chroms, "start":pos, "variant_id":variant_ids, "AF":AFs, "end":ends, "type":types}
df = pd.DataFrame(dict_)
df = df.reindex(columns=['chr', 'start', 'end','type', 'variant_id', 'AF'])
df.to_csv(outpath, index=False, sep="\t", header=False)
os.system("sort -V -k1,1 -k2,2 "+outpath+" > "+os.path.splitext(outpath)[0]+".sorted.tsv")


