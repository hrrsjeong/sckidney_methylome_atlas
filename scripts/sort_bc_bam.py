import sys,glob
import os
#filtered_bc_bam/kidney-atlas_ACCTCAT/kidney-atlas_TCTATCGGTA+ACTATAGGTT+AGCCTCAT.bam
input_dir_flag = sys.argv[1] #"filtered_bc_bam/kidney-atlas_ACCTCAT/split.complete"#sys.argv[1]
if not os.path.exists(input_dir_flag):
    raise ValueError("bam split by bc not completed correctly")
input_dir = '/'.join(input_dir_flag.split('/')[0:-1])
complete_flag = sys.argv[2] # "filtered_bc_bam/kidney-atlas_ACCTCAT/sorted/sort.complete"#sys.argv[2]
output_dir = '/'.join(complete_flag.split('/')[0:-1])
if os.path.exists(output_dir):
    os.system(f"rm -rf {output_dir}") 
os.system(f"mkdir -p {output_dir}")
flist = glob.glob(f'{input_dir}/*.bam')
for sfile in flist:
    sfile_out_base = sfile.split("/")[-1].replace(".bam",".sorted.bam")
    systemstr = f'samtools sort -@2 -o {output_dir}/{sfile_out_base} {sfile} && samtools index {output_dir}/{sfile_out_base}'
    #print (systemstr)
    os.system(systemstr)
os.system(f"touch {complete_flag}")
