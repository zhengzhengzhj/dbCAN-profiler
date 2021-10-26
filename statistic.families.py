### 
import sys

# python3 statistic.families.py HMP_rundbcan/data/dbcan_out/overview.txt

def output(hmmcount):
    for p in hmmcount:
        print (f"{p}\t{hmmcount[p]}")

def HMM_f(filename):
    hmmcount = {}
    for line in open(filename).readlines()[1:]:
        lines = line.split()
        pfams = lines[1].split("+")
        for i in pfams:
            a = i.split("(")[0]
            a = a.split("_")[0]
            if a == "-":
                continue
            if a in hmmcount:
                hmmcount[a] += 1
            else:
                hmmcount[a] = 1
    output(hmmcount)   

HMM_f(sys.argv[1])

