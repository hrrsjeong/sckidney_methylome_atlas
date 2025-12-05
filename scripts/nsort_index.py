import sys,glob,os
import pysam
os.system("mkdir -p tmp2/")
fin = "rmdup_bam/kidney-atlas.AGCCTCAT.bbrd.q10.bam" #sys.argv[1]
ftmp = f"tmp2/{fin.split('/')[-1]}"
ftmp = ftmp.replace(".bam",".nsort.bam")
os.system(f"samtools sort -n -@4 -o {ftmp} {fin}")
#os.system(f"samtools index {ftmp}")

samfile = pysam.AlignmentFile(ftmp, "rb")
#header = samfile.header
#samfile.close()


