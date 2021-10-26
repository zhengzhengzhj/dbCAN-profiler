import os
command = 'seqkit stats sample0_P1/data/input1_val_1.fq.gz'#可以直接在命令行中执行的命令
r = os.popen(command)
info = r.readlines()
count = 0
for line in info:
    line = line.strip('\n')
    count += 1
    print (count,line)
