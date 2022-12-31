import sys
from collections import Counter
cnt = Counter()
with open(sys.argv[1],'r') as f:
    for line in f:
        line = line.rstrip().split()
        gene = line[0] + '|' + line[2]
        cnt[gene]+=1

for k,v in cnt.items():
    print(k,v)
