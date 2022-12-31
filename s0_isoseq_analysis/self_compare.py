import sys
from cupcake.tofu.compare_junctions import compare_junctions
from cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format
from collections import defaultdict

PB_gene = {}
with open(sys.argv[1],'r') as f:
    for line in f:
        PB, gene = [line.rstrip().split()[i] for i in [0,11]]
        PB_gene[PB] = gene


recs = defaultdict(lambda: [])
ref = defaultdict(lambda: [])
recs_reader = [r for r in collapseGFFReader(sys.argv[2])]
for r in recs_reader:
    recs[PB_gene[r.seqid]].append(r)

ref_reader =  [r for r in collapseGFFReader(sys.argv[2])]
for r in ref_reader:
    ref[PB_gene[r.seqid]].append(r)

for g in recs:
    for i in recs[g]:
        for j in ref[g]:
            if i.seqid != j.seqid:
                i.segments = i.ref_exons
                j.segments = j.ref_exons
                match_type = compare_junctions(i, j, internal_fuzzy_max_dist=10, max_5_diff=5000, max_3_diff=50)
                print (str(g)+'\t'+i.seqid+'\t'+j.seqid+'\t'+match_type)
