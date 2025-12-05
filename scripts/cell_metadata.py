import sys,glob,os

#Tn5_list = sys.argv[1] #Tn5_bismark_extract_dir.kidney-atlas.txt
cell_annotation = sys.argv[2] #"/mnt/nfs/home/Projects/sciMet/NovaSeqS4_230501/scaleMethyl_outs/Sauron.annot"
cell_metadata = sys.argv[3] #outfile
#filtered_bc_bismark2extract/kidney-atlas_AACGCATC/bismark_extract.complete

allcools_file_list = sys.argv[1] #allcools_file_list.kidney-atlas.txt
#allcools_file_list.kidney-atlas.txt
'''
head filtered_bc_CXreport/kidney-atlas_ATGAGGCC/kidney-atlas_GCCATCAACT+TGCTAATTCT+ATGAGGCC.cytosine_context_summary.txt
upstream        C-context       full context    count methylated        count unmethylated      percent methylation
A       CAA     ACAA    7862    1064067 0.73
C       CAA     CCAA    7323    793390  0.91
G       CAA     GCAA    5143    713058  0.72
T       CAA     TCAA    6341    984797  0.64

head allcools_file_list.kidney-atlas2022.txt 
filtered_bc_allcools/kidney-atlas2022_AACGTTCC/kidney-atlas2022_AGTTGATGGC+GGTACCGGCT+AACGTTCC.deduplicated.allc.tsv.gz
filtered_bc_allcools/kidney-atlas2022_AACGTTCC/kidney-atlas2022_TCCAGAAGCC+TTGGTCGAAT+AACGTTCC.deduplicated.allc.tsv.gz

head filtered_bc_CXreport/kidney-atlas2022_ATAATGTG/kidney-atlas2022_GAAGAGTATT+CCAATTCCAT+ATAATGTG.deduplicated.cytosine_context_summary.txt
upstream	C-context	full context	count methylated	count unmethylated	percent methylation
A	CAA	ACAA	1176	208949	0.56
C	CAA	CCAA	1526	177301	0.85
G	CAA	GCAA	994	155287	0.64
T	CAA	TCAA	1439	206922	0.69


GCAGGTCGTCACGAGGAGCCTATTAGCT	K19
GCAGGTCGTCACGAGGAGCCGGAGCGTC	K19
head filtered_bc_allcools/kidney-atlas2022_AACGTTCC/kidney-atlas2022_AGTTGATGGC+GGTACCGGCT+AACGTTCC.deduplicated.allc.tsv.gz

cat filtered_bc_rmdup_bam/kidney-atlas2022_AGTACTCC/kidney-atlas2022_ACGCGACGGC+GTAGTAGTCC+AGTACTCC.deduplication_report.txt

Total number of alignments analysed in filtered_bc_bam/kidney-atlas2022_AGTACTCC/kidney-atlas2022_ACGCGACGGC+GTAGTAGTCC+AGTACTCC.bam:	99189
Total number duplicated alignments removed:	5720 (5.77%)
Duplicated alignments were found at:	4575 different position(s)

Total count of deduplicated leftover sequences: 93469 (94.23% of total)

'''
dd = {}
rmdup = set()
with open(allcools_file_list,'r') as fp:
    for line in fp:
        sample_Tn5 = line.strip().split('/')[-2]
        sample,Tn5 = sample_Tn5.split('_')
        if sample_Tn5 in rmdup:continue
        rmdup.add(sample_Tn5)
        cytosine_count = glob.glob(f"filtered_bc_CXreport/{sample_Tn5}/*.cytosine_context_summary.txt")#count.csv")
        for cell_count_file in cytosine_count:
            cell_name = cell_count_file.split('/')[-1].split('.')[0]
            dd.setdefault(cell_name,{})
            dd[cell_name].setdefault("c_count",{"CH":[0,0],"CG":[0,0],"CCC":[0,0]})
            dd[cell_name].setdefault("c_frac",{"CH":0,"CG":0,"CCC":0})
            with open(cell_count_file,'r') as fcell:
                fcell.readline()
                for line in fcell:
                    line_temp = line.strip().split('\t')
                    c_base = line_temp[1]
                    mc,cov,mc_rate = int(line_temp[3]),int(line_temp[3])+int(line_temp[4]),float(line_temp[5])/100
                    if "N" in c_base:continue
                    if c_base.startswith("CG"):
                        dd[cell_name]["c_count"]["CG"][0] += mc
                        dd[cell_name]["c_count"]["CG"][1] += cov
                    else:
                        dd[cell_name]["c_count"]["CH"][0] += mc
                        dd[cell_name]["c_count"]["CH"][1] += cov
                    if c_base == "CCC":
                        dd[cell_name]["c_count"]["CCC"][0] += mc
                        dd[cell_name]["c_count"]["CCC"][1] += cov
                        #dd[cell_name]["c_frac"]["CCC"] = mc_rate
            dd[cell_name]["c_frac"]["CG"] = float(dd[cell_name]["c_count"]["CG"][0]) / float(dd[cell_name]["c_count"]["CG"][1])
            dd[cell_name]["c_frac"]["CH"] = float(dd[cell_name]["c_count"]["CH"][0]) / float(dd[cell_name]["c_count"]["CH"][1])
            dd[cell_name]["c_frac"]["CCC"] = float(dd[cell_name]["c_count"]["CCC"][0]) / float(dd[cell_name]["c_count"]["CCC"][1])
        with open(f"temp_complexity/{sample}.{Tn5}.complexity.txt",'r') as comp:
            for line in comp:
                line_temp = line.strip().split('\t')
                cell_id = f"{sample}_{line_temp[1]}"
                #total_read,uniq_read,uniq_read_percent = line_temp[2:]
                if cell_id in dd:
                    dd[cell_id]["read"] = line_temp[2:]
        rmdup_list_per_Tn5 = glob.glob(f"filtered_bc_rmdup_bam/{sample}_{Tn5}/{sample}_*.deduplication_report.txt")
        for rmdup_cell_file in rmdup_list_per_Tn5:
            cell_id = rmdup_cell_file.split('/')[-1].split('.dedupli')[0]
            with open(rmdup_cell_file,'r') as fin_rmdup_cell_file:
                flag = 0
                for line in fin_rmdup_cell_file:
                    #line_temp = line.strip().split('\t')
                    #cell_id = f"{sample}_{line_temp[1]}"
                    #total_read,uniq_read,uniq_read_percent = line_temp[2:]
                    if line.startswith("Total count of"):
                        flag = 1
                        tcount = line.strip().split(": ")[-1].split(' ')[0]
                        dd[cell_id]["rmdup_read"] = tcount
                if flag == 0:
                    print ("dedup report not available: ",rmdup_cell_file)
'''
Total number of alignments analysed in filtered_bc_bam/kidney-atlas_AGACTCGA/kidney-atlas_CGAATCTCCT+CAACCAGTAC+AGACTCGA.bam:   1041909
Total number duplicated alignments removed:     109327 (10.49%)
Duplicated alignments were found at:    78973 different position(s)

Total count of deduplicated leftover sequences: 932582 (89.51% of total)
'''
                
with open(cell_annotation,'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        cell_idx1,cell_idx2,cell_idx3 = line_temp[0][0:10],line_temp[0][10:20],line_temp[0][20:]
        cell_id = f"{sample}_{cell_idx1}+{cell_idx2}+{cell_idx3}"
        if cell_id in dd:
            dd[cell_id]["sample_id"] = line_temp[1]

with open(cell_metadata,'w') as fout:
    fout.write("index\tmCCCFrac\tmCGFrac\tmCHFrac\tTotalReads\tUniqReads_tmp\tUniqReadsPercent\tUniqReads\tCG_overlap_cnt\tSampleID\n")
    for cell_id in dd:
        if "rmdup_read" not in dd[cell_id]:
            print (cell_id)
            continue
        fout.write(f'{cell_id}\t{dd[cell_id]["c_frac"]["CCC"]}\t{dd[cell_id]["c_frac"]["CG"]}\t{dd[cell_id]["c_frac"]["CH"]}')
        fout.write('\t'+'\t'.join(dd[cell_id]["read"]))
        fout.write('\t'+dd[cell_id]["rmdup_read"])
        fout.write(f'\t{dd[cell_id]["c_count"]["CG"][1]}')
        fout.write('\t'+dd[cell_id]["sample_id"])
        fout.write('\n')


'''
TTGATATGAACTTCTCATTGTATCAATC    kidney-b
TTGATATGAACTTCTCATTGGCATGAAT    kidney-b
TTGATATGAACTTCTCATTGAACTTCCA    kidney-b
TTGATATGAA CTTCTCATTG TTGGTACC    kidney-b
AGCGATCCGC+TTGCCTTGGC+AGGTGATT  10286934        4893287 47.57
~/data/net/home/Projects/sciMet/processing/Human_kidney_atlas_Novaseq  more ../../NovaSeqS4_230501/scaleMethyl_outs/Sauron.annot 
'''

                    
                
        





'''
head complexity/kidney-atlas.AGGTGATT.complexity.txt 
1       AGCGATCCGC+TTGCCTTGGC+AGGTGATT  10286934        4893287 47.57
2       ATATGCCATC+TAACGAATTG+AGGTGATT  12779344        4830880 37.80
3       AGCGATCCGC+TGAGAACCAA+AGGTGATT  11047332        4396944 39.80

filtered_bc_allc/kidney-atlas_ACCTTGGC/kidney-atlas_AAGCATCCTA+CATAGTTCGG+ACCTTGGC.sorted.gz.count.csv
,mc,cov,mc_rate,genome_cov
CAC,3949,580076,0.006807728642453747,0.012349521133331341
CAG,4290,810576,0.0052925327174749804,0.012349521133331341
CTC,1689,659734,0.0025601227161249837,0.012349521133331341
CGT,74100,95566,0.7753803654019212,0.012349521133331341
CTA,1408,558372,0.002521616413430473,0.012349521133331341
CAT,4314,797932,0.00540647573978735,0.012349521133331341
CGC,60725,76363,0.7952149601246677,0.012349521133331341
CAA,3840,836467,0.0045907369926129786,0.012349521133331341
CGA,72006,91303,0.7886487848153949,0.012349521133331341
CCG,415,92853,0.0044694301745770194,0.012349521133331341
CCA,2052,742151,0.0027649359766408724,0.012349521133331341
CCC,983,481673,0.0020408036157310043,0.012349521133331341
CTG,2162,818106,0.00264268933365603,0.012349521133331341
CCT,1700,703377,0.002416911556675865,0.012349521133331341
CTT,2190,870607,0.0025154863216123925,0.012349521133331341
CGG,70707,91981,0.7687131037931747,0.012349521133331341
CTN,0,1,0.0,0.012349521133331341
CNN,0,3,0.0,0.012349521133331341
'''
