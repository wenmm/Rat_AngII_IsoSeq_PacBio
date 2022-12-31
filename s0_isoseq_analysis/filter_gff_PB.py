import sys
import re

delete_PB = []
with open(sys.argv[1],'r') as f:
    for line in f:
        pb = line.rstrip().split()[0]
        delete_PB.append(pb)

regexTID  = re.compile ('transcript_id \"([^\"]+)\"\;')
with open(sys.argv[2],'r') as f:
    for line in f:
        lines = line.strip().split("\t")
        attrs = lines[8]
        transcriptID = re.search(regexTID, attrs).group(1)
        if transcriptID in delete_PB:
            pass
        else:
            print (line.rstrip())
