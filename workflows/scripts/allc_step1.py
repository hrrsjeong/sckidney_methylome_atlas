import sys,glob
import os
import queue,threading

class ThreadRBH(threading.Thread):
    def __init__(self,queue):
        threading.Thread.__init__(self)
        self.queue = queue
    def run(self):
        while True:
            qcommand = self.queue.get()
            os.system(qcommand)
            self.queue.task_done()

exe_list = []

#filtered_bc_bam/kidney-atlas_ACCTCAT/kidney-atlas_TCTATCGGTA+ACTATAGGTT+AGCCTCAT.bam
ref_file_path = sys.argv[1]
input_dir_flag = sys.argv[2] #"filtered_bc_bam/sorted/{sample}_{Tn5}/sort.complete"#sys.argv[1]
if not os.path.exists(input_dir_flag):
    raise ValueError("sorted bc bam not completed correctly")
input_dir = '/'.join(input_dir_flag.split('/')[0:-1])
complete_flag = sys.argv[3] # "filtered_bc_allc/{sample}_{Tn5}/allc_step1.complete"#sys.argv[2]
mapq = sys.argv[4]
output_dir = '/'.join(complete_flag.split('/')[0:-1])
if os.path.exists(output_dir):
    os.system(f"rm -rf {output_dir}")
os.system(f"mkdir -p {output_dir}")
flist = glob.glob(f'{input_dir}/*.bam')
for sfile in flist:
    sfile_out_base = sfile.split("/")[-1].split(".bam")[0]
    systemstr = f'allcools bam-to-allc --convert_bam_strandness --bam_path {sfile} --save_count_df --min_mapq {mapq} --reference_fasta {ref_file_path} --output_path {output_dir}/{sfile_out_base}.allc.tsv'
    exe_list.append(systemstr)
    #print (systemstr)
    #os.system(systemstr)

queue = queue.Queue()
for i in range(4):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
    for exes in exe_list:
        queue.put(exes)
        queue.join()

os.system(f"touch {complete_flag}")
