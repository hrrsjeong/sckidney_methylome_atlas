import sys
fout = open(sys.argv[2],'w')
with open(sys.argv[1],'r') as fp:
    for line in fp:
        line_temp = line.strip().split('/')[-1].split('.')[0]
        fout.write(line_temp+'\t'+line.strip()+'\n')
fout.close()
