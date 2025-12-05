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

ref_file_path = sys.argv[1]
ref_folder = '/'.join(ref_file_path.split('/')[0:-1])
input_dir_flag = sys.argv[2] #"filtered_bc_bismark2extract/{sample}_{Tn5}/bismark_extract.complete"#sys.argv[1]
if not os.path.exists(input_dir_flag):
    raise ValueError("sorted bc bam not completed correctly")
input_dir = '/'.join(input_dir_flag.split('/')[0:-1])
complete_flag = sys.argv[3] # "filtered_bc_CXreport/{sample}_{Tn5}/coverage2cytosine.complete"#sys.argv[2]
output_dir = '/'.join(complete_flag.split('/')[0:-1])
if os.path.exists(output_dir):
    os.system(f"rm -rf {output_dir}")
os.system(f"mkdir -p {output_dir}")
flist = glob.glob(f'{input_dir}/*.bismark.cov.gz')
for sfile in flist:
    sfile_out_base = sfile.split("/")[-1].split(".bismark.cov.gz")[0] #../filtered_bc_bismark2extract/kidney-atlas_GGCATTCT/kidney-atlas_ACGCTTCTCT+TGCTAATTCT+GGCATTCT.bismark.cov.gz
    systemstr = f'coverage2cytosine --coverage_threshold 1 --CX_context --genome_folder {ref_folder} --dir {output_dir} --gzip -o {sfile_out_base} {sfile}'
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
