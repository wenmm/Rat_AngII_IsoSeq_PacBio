import sys
from bisect import bisect

ref_pas = {}
polya_db = {}

with open(sys.argv[1],'r') as f:
    for line in f:
        raw = line.rstrip().split()
        ref_pas.setdefault(raw[0],[]).append(int(raw[1]))

with open(sys.argv[2],'r') as f:
    for line in f:
        raw = line.rstrip().split()
        polya_db.setdefault(raw[0]+'|'+raw[5],[]).append(int(raw[1]))

with open(sys.argv[3],'r') as f:
#    next(f)
    for line in f:
        raw = line.rstrip().split()
        tag = raw[1] + '|' + raw[2]
        pos = int(raw[3])
        strand = raw[2]
        try:
            a = bisect(ref_pas[tag],pos)
            b = bisect(polya_db[tag],pos)
            #print(ref_pas[tag][a])
            print(a,b,ref_pas[tag][a],b,polya_db[tag])
            if a == len(ref_pas[tag]):
                near_ref = ref_pas[tag][a-1]
            else:
                near_ref = ref_pas[tag][a]
            if b == len(polya_db[tag]):
                near_polya_db = polya_db[tag][b-1]
            else:
                near_polya_db = polya_db[tag][b]
            near_ref_diff = abs(ref_pas[tag][a] - pos)
            near_polya_db_diff = abs(polya_db[tag][b] - pos)
            if abs(ref_pas[tag][a] - pos) >= abs(ref_pas[tag][a-1] - pos):
                near_ref = ref_pas[tag][a-1]
            if abs(polya_db[tag][b] - pos) >= abs(polya_db[tag][b-1] - pos):
                near_polya_db = polya_db[tag][b-1]
            diff1 = abs(near_ref - pos)
            diff2 = abs(near_polya_db - pos)
            if diff1 <= diff2:
                queding_pos = near_ref
            else:
                queding_pos = near_polya_db
            if strand == '+':
                print("{}\t{}\t{}".format(abs(pos - queding_pos) , queding_pos, line.rstrip()))
            if strand == '-':
                print("{}\t{}\t{}".format(abs(queding_pos - pos),queding_pos, line.rstrip()))

        except:
            pass
#        print("{}\t{}\t{}",format(a, b, line.rstrip()))
