import sys,os,glob
import queue,threading

class ThreadRBH(threading.Thread):
    def __init__(self,queue):
        threading.Thread.__init__(self)
        self.queue = queue
    def run(self):
        while True:
            systemstr = self.queue.get()
            print (systemstr)
            os.system(systemstr)
            self.queue.task_done()

exe_list = []
#sample_folder_tmp = [x for x in sys.argv[1].split('/') if x.strip() != ""]
#python /mnt/nfs/home/Projects/sciMet/scimetv2_pipeline/scripts/sort_chroms.py kidney-A.CG.chrom_folders.txt 16 kidney-A.CG.chrom_folders.sort.done
sample_folder_path = sys.argv[1]

flist = []
with open(sample_folder_path,'r') as fp:
    for line in fp:
        sample_folder = line.strip()
        print (sample_folder)
        flist_tmp = glob.glob(f"{sample_folder}/*.bed")
        flist += flist_tmp
print (sample_folder_path)
print (flist)

for sfile in flist:
    tmp_systemstr = f"cat {sfile} | sort -k2,2n | gzip >| {sfile}.gz"
    exe_list.append(tmp_systemstr)

queue = queue.Queue()
for i in range(int(int(sys.argv[2])/2)):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
queue.join()
#sort_flag = sample_folder_path.replace('.txt','.sort.done')
os.system(f"echo done >| {sys.argv[3]}")

'''
perl /mnt/nfs/home/Projects/sciMet/scimetv2_pipeline/sciMETv2/sciMET_sortChroms.pl kidney-A.CG.011.chroms
Command:
	for f in kidney-A.CG.011.chroms/*.bed; do cat $f | sort -k2,2n | gzip > $f.gz & done
Copy/paste then run command.

  hjeong@ip-10-112-107-181  ❲c❳ test  ~/data/net/home/Projects/sciMet/processing_tmp  ls kidne
Display all 127 possibilities? (y or n)
  hjeong@ip-10-112-107-181  ❲c❳ test  ~/data/net/home/Projects/sciMet/processing_tmp  ls kidney-A.CG.011.chroms/
chr1.bed   chr11.bed  chr13.bed  chr15.bed  chr17.bed  chr19.bed  chr20.bed  chr22.bed  chr4.bed  chr6.bed  chr8.bed  chrX.bed
'''
