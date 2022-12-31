import sys
from cupcake.tofu.compare_junctions import compare_junctions
from cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format
from collections import defaultdict

gene_id = {}
with open(sys.argv[1],'r') as f:
    for line in f:
        raw = line.rstrip().split()
        gene_id.setdefault(raw[0],[]).append(raw[1])

recs = defaultdict(lambda: [])
recs_reader = [r for r in collapseGFFReader(sys.argv[2])]

ref = defaultdict(lambda: [])
ref_reader = [r for r in collapseGFFReader(sys.argv[3])]

for r in recs_reader:
    if r.seqid in gene_id:
        for h in gene_id[r.seqid]:
            recs[h].append(r)
    
for r in ref_reader:
    ref[r.geneid].append(r)

for g in recs:
    for i in recs[g]:
        for j in ref[g]:
            if i.seqid != j.seqid:
                i.segments = i.ref_exons
                j.segments = j.ref_exons
                match_type = compare_junctions(i, j, internal_fuzzy_max_dist=10, max_5_diff=50000000, max_3_diff=50)
                print (str(g)+'\t'+i.seqid+'\t'+j.seqid+'\t'+match_type)
            
