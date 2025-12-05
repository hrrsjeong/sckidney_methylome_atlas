import sys
fout_pass_cell = open(sys.argv[2],'w')
num_uniq_bc = int(sys.argv[3])
with open(sys.argv[1],'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        cell_id = line_temp[1]
        if int(line_temp[3]) > num_uniq_bc:
            fout_pass_cell.write(cell_id+'\n')
fout_pass_cell.close()




'''
cat kidney-A.merged.complexity.txt | head
1       AATTGAGAGA+GTAGCTCCAT+GACTTACA  9184914 5283736 57.53
2       AATTGAGAGA+TGCTAATTCT+AACTGTAG  6706276 3920198 58.46

'''
