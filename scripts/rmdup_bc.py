import sys,glob
import os
#filtered_bc_bam/kidney-atlas_ACCTCAT/kidney-atlas_TCTATCGGTA+ACTATAGGTT+AGCCTCAT.bam
input_dir_flag = sys.argv[1] #"filtered_bc_bam/kidney-atlas_ACCTCAT/split.complete"#sys.argv[1]
input_dir = '/'.join(input_dir_flag.split('/')[0:-1])
complete_flag = sys.argv[2] # "filtered_bc_bismark2extract/kidney-atlas_ACCTCAT/bismark_extract.complete"#sys.argv[2]
output_dir = '/'.join(complete_flag.split('/')[0:-1])
os.system(f"mkdir -p {output_dir}")
flist = glob.glob(f'{input_dir}/*.bam')
for sfile in flist:
    systemstr = f'deduplicate_bismark -p --output_dir {output_dir} {sfile}' #filtered_bc_bam/kidney-atlas_ATGAGGCC/*.bam
    print (systemstr)
    os.system(systemstr)
#systemstr = f'deduplicate_bismark -p --output_dir {output_dir} {input_dir}/*.bam' #filtered_bc_bam/kidney-atlas_ATGAGGCC/*.bam
#print (systemstr)
#os.system(systemstr)
os.system(f"touch {complete_flag}")
