import sys,glob
import os
import gzip

#filtered_bc_bam/kidney-atlas_ACCTCAT/kidney-atlas_TCTATCGGTA+ACTATAGGTT+AGCCTCAT.bam
input_dir_flag = sys.argv[1] #"filtered_bc_CXreport/{sample}_{Tn5}/CXreport.complete"#sys.argv[1]
'''
#GL000008.2      1078    +       0       1       CHH     CAC
GL000008.2      1089    +       0       2       CHH     CCT
GL000008.2      1090    +       0       2       CHG     CTG
GL000008.2      1094    +       1       1       CHH     CAA
 /mnt/nfs/home/Projects/sciMet/processing/Human_kidney_atlas_Novaseq_3/filtered_bc_CXreport/kidney-atlas_TACTGGTC/kidney-atlas_CCTAAGCGGT+CATAGTTCGG+TACTGGTC.CX_report.txt.gz
''' 
if not os.path.exists(input_dir_flag):
    raise ValueError("CX report inputs not completed correctly")

input_dir = '/'.join(input_dir_flag.split('/')[0:-1])
#complete_flag = sys.argv[3] # "filtered_bc_allc/{sample}_{Tn5}/filtered/allc_step2.complete"#sys.argv[2]
output_count = sys.argv[2]
output_dir = '/'.join(output_count.split('/')[0:-1])
#if os.path.exists(output_dir):
#    os.system(f"rm -rf {output_dir}")
#os.system(f"mkdir -p {output_dir}/")
flist = glob.glob(f'{input_dir}/*.CX_report.txt.gz')
fout = open(output_count,'w')
for sfile in flist:
    sfile_out_base = sfile.split("/")[-1].split(".CX_report.txt.gz")[0]
    dd = {}
    with gzip.open(sfile,'rt') as fp:
        for line in fp:
            line_temp = line.strip().split('\t')
            c_base = line_temp[5]
            cnt = int(line_temp[3])+int(line_temp[4])
            dd.setdefault(c_base,{})
            dd[c_base].setdefault(cnt,0)
            dd[c_base][cnt] +=1
    for c_base in dd:
        for cnt in dd[c_base]:
            fout.write(f"{sfile_out_base}\t{c_base}\t{cnt}\t{dd[c_base][cnt]}\n")
fout.close()


'''
CXreport
GL000008.2      474     +       0       3       CHH     CCT
GL000008.2      475     +       0       1       CHH     CTT

allc
chr1    52005   +       CTG     0       1       1



python /mnt/nfs/home/Program/pipelines/scimetv2_pipeline/scripts/allc_step2_rmSNP.py /mnt/nfs/home/Projects/sciMet/ref/dbSnp155Common.snp.bed filtered_bc_allc/kidney-atlas_AATGCCTC/allc_step1.complete filtered_bc_allc/autosome_rmSNP/kidney-atlas_AATGCCTC/allc_step2.complete
queue = queue.Queue()
for i in range(4):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
    queue.join()
'''
#os.system(f"rm -rf {output_dir}/tmp_allc_sh")
#os.system(f"touch {complete_flag}")
