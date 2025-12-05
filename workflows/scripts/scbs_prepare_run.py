import sys,os
import glob

input_dir_flag = sys.argv[1]
output_dir_flag = sys.argv[2]

output_dir = '/'.join(output_dir_flag.split('/')[0:-1])
cov_dir_list = []
with open(input_dir_flag,'r') as fp:
    for line in fp:
        tmp_dir_path = '/'.join(line.strip().split('/')[0:-1])+'/*.cov.gz'
        tmp_dir_check = glob.glob(tmp_dir_path)
        if len(tmp_dir_check) < 1:continue
        cov_dir_list.append(tmp_dir_path)
systemstr = "scbs prepare "+" ".join(cov_dir_list)+" "+output_dir
print (systemstr)
os.system(systemstr)

