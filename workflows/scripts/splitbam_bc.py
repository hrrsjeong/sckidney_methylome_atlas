import sys,glob,os
import pysam
fin = sys.argv[1] #"rmdup_bam/kidney-atlas.AGCCTCAT.bbrd.q10.bam" #sys.argv[1]
bc_file = sys.argv[2] #"complexity/kidney-atlas.AGCCTCAT.complexity.txt" #sys.argv[2]
complete_flag = sys.argv[3] #"filtered_bc_bam/kidney-atlas_ACCTCAT/kidney-atlas" #sys.argv[3]
out_dir = '/'.join(complete_flag.split('/')[0:-1])
base_out = sys.argv[4]#'/'.join(base_out.split('/')[0:-1])
os.system(f"mkdir -p {out_dir}")
num_uniq_reads = int(sys.argv[5]) #200000
filtered_bc_files = {}
samfile = pysam.AlignmentFile(fin, "rb",require_index=False)
header = samfile.header

with open(bc_file,'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        if int(line_temp[3]) > num_uniq_reads:
            filtered_bc = line_temp[1]
            filtered_bc_filename = f"{out_dir}/{base_out}_{filtered_bc}.bam"
            bam_out = pysam.AlignmentFile(filtered_bc_filename, "wb", header=header)
            filtered_bc_files[filtered_bc] = bam_out
flag = 0
pre_readID = ''
pre_info = ''
for align in samfile:
    readID = align.query_name
    bcID = readID.split(':')[0]
    if bcID not in filtered_bc_files:continue
    if pre_readID == readID:
        filtered_bc_files[bcID].write(pre_info)
        filtered_bc_files[bcID].write(align)
    pre_readID = readID
    pre_info   = align

samfile.close()

for bc in filtered_bc_files:
    filtered_bc_files[bc].close()

os.system(f"touch {complete_flag}")
'''
1	GGCTTCTGGA+CAGGAGGAGA+AGCCTCAT	1553130	746915	48.09
2	TGGCAGAAGT+CATCCGACTA+AGCCTCAT	1535316	659524	42.96
3	TCTATCGGTA+ACTATAGGTT+AGCCTCAT	1186038	634470	53.49
'''
