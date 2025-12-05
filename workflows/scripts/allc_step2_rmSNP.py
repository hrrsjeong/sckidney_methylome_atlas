import sys,glob
import os
#import queue,threading
'''
class ThreadRBH(threading.Thread):
    def __init__(self,queue):
        threading.Thread.__init__(self)
        self.queue = queue
    def run(self):
        while True:
            qcommand = self.queue.get()
            os.system(qcommand)
            self.queue.task_done()
'''
exe_list = []

#filtered_bc_bam/kidney-atlas_ACCTCAT/kidney-atlas_TCTATCGGTA+ACTATAGGTT+AGCCTCAT.bam
SNP_file_path = sys.argv[1]
input_dir_flag = sys.argv[2] #"filtered_bc_allc/{sample}_{Tn5}/allc_step1.complete"#sys.argv[1]
if not os.path.exists(input_dir_flag):
    raise ValueError("allc step1 input not completed correctly")

input_dir = '/'.join(input_dir_flag.split('/')[0:-1])
complete_flag = sys.argv[3] # "filtered_bc_allc/{sample}_{Tn5}/filtered/allc_step2.complete"#sys.argv[2]
output_dir = '/'.join(complete_flag.split('/')[0:-1])
if os.path.exists(output_dir):
    os.system(f"rm -rf {output_dir}")
os.system(f"mkdir -p {output_dir}/tmp_allc_sh/")
flist = glob.glob(f'{input_dir}/*.allc.tsv.gz')
for sfile in flist:
    sfile_out_base = sfile.split("/")[-1].split(".allc.tsv.gz")[0]
    fout_sh_name = f'{output_dir}/tmp_allc_sh/{sfile_out_base}.sh'
    with open(fout_sh_name,'w') as fout_sh:
        fout_sh.write("zcat "+sfile+" | awk -v OFS='\\t' '{print $1,$2-1,$2,$3,$4,$5,$6,$7}' | bedtools subtract -a - -b "+SNP_file_path+' | cut -f 1,3,4,5,6,7,8 | grep -Fv -e "chrX" -e "chrY" -e "chrM" -e "KI2" -e "GL0"  | bgzip -c > '+output_dir+'/'+sfile_out_base+'.allc.tsv.gz && tabix -b 2 -e 2 -s 1 '+output_dir+'/'+sfile_out_base+'.allc.tsv.gz\n')
    systemstr = f"bash {fout_sh_name}"
    os.system(systemstr)
    #exe_list.append(systemstr)
'''
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
os.system(f"rm -rf {output_dir}/tmp_allc_sh")
os.system(f"touch {complete_flag}")
