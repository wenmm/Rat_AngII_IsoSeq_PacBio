import sys
import math

control_gene_count = {}
angii_gene_count = {}

with open(sys.argv[1],'r') as f:
    for line in f:
        raw = line.rstrip().split()
        control_gene_count.setdefault(raw[0],[]).append(int(raw[2]))
        angii_gene_count.setdefault(raw[0],[]).append(int(raw[3]))

with open(sys.argv[1],'r') as f:
    for line in f:
        raw = line.rstrip().split()
        if sum(control_gene_count[raw[0]]) > 0 and sum(angii_gene_count[raw[0]]) > 0:
            fc = math.log2(sum(angii_gene_count[raw[0]])/sum(control_gene_count[raw[0]]))
            print(line.rstrip() + '\t' + str(int(raw[2])/sum(control_gene_count[raw[0]])) + '\t' + str(int(raw[3])/sum(angii_gene_count[raw[0]])) + '\t' + str(fc))
        if sum(control_gene_count[raw[0]]) == 0 and sum(angii_gene_count[raw[0]]) > 0:
            print(line.rstrip() + '\t' + '0' +  '\t' + str(int(raw[3])/sum(angii_gene_count[raw[0]])) + '\t' + '999999999999')
        if sum(control_gene_count[raw[0]]) > 0 and sum(angii_gene_count[raw[0]]) == 0:
            print(line.rstrip() + '\t' + str(int(raw[2])/sum(control_gene_count[raw[0]])) + '\t' + '0' + '\t' + '-999999999999')

