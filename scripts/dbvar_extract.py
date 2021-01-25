import gzip
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-vcf', '--vcf',
                    help='Path to the dbvar_SV vcf file', required=True)
parser.add_argument('-chr', '--chr',
                    help='the chromosome', required=True)
parser.add_argument('-out', '--out',
                    help='Path to the output tsv file', required=True)
args = parser.parse_args()

outf = open(args.out, 'w')
# header
outf.write('chr\tstart\tend\ttype\tvariant_id\n')
# read VCF with all variants
for line in gzip.open(args.vcf, 'rb'):
    line = line.decode('ascii')
    # skip headers
    if line[0] == '#':
        continue
    line = line.rstrip().split('\t')
    # select chromosome
    if line[0] != args.chr:
        continue
    # parse INFO field
    infos = {}
    for info in line[7].split(';'):
        info_pair = info.split('=')
        if len(info_pair) == 1:
            infos[info] = True
        else:
            infos[info_pair[0]] = info_pair[1]
    # prepare info for TSV
    out_line = [line[0]]
    # start pos
    pos_s = line[1]
    out_line.append(pos_s)
    # end pos
    pos_e = pos_s
    if 'END' in infos:
        pos_e = infos['END']
    out_line.append(pos_e)
    # SV type
    svtype = infos['SVTYPE']
    out_line.append(svtype)
    # skip if BND type
    if svtype == 'BND':
        continue
    if int(pos_e) < int(pos_s):
        print(svtype)
    # variant id
    out_line.append(line[2])
    # write line in tsv
    outf.write('\t'.join(out_line) + '\n')
outf.close()
