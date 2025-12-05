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
    systemstr = f'bismark_methylation_extractor --comprehensive --merge_non_CpG --CX_context --no_header --gzip --multicore 2 --no_overlap --bedGraph --buffer_size 24G -o {output_dir} {sfile}'
    print (systemstr)
    os.system(systemstr)
os.system(f"touch {complete_flag}")
