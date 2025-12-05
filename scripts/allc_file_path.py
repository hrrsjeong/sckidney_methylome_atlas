import sys,glob

with open(sys.argv[2],'w') as fout:
    with open(sys.argv[1],'r') as fp: #input is complete flag
        for line in fp:
            tmp_match = '/'.join(line.strip().split('/')[0:-1])+'/*.allc.tsv.gz'
            tmp_flist = glob.glob(tmp_match)
            for i in tmp_flist:
                fout.write(i+'\n')
