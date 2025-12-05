import pandas as pd
#import Bio.SeqIO
#import Bio.Seq
#import collections
import gzip
import numpy as np
import os
#import pysam
import re
import shutil

SDIR=os.path.dirname(workflow.snakefile)
CWD=os.getcwd()
#shell.prefix(f"source {SDIR}/env.cfg ; set -eo pipefail; ")
#shell.prefix(f"source ~/.bashrc")
configfile: 'config.yaml'

min_map_score = config['MIN_MAP_SCORE']
num_uniq_reads = config['NUM_UNIQ_READS']
#sliding_win_file = config['SLIDING_WIN_FILE']
reference_fasta = config["REF"]
bsbolt_dir = config["REF"].split('/')[-1].split('.fa')[0].split('.fasta')[0]
reference_dir = '/'.join(reference_fasta.split('/')[0:-1])
reference_fasta_fai = reference_fasta+'.fai'
if not os.path.exists(reference_fasta_fai):
	os.system(f"samtools faidx {reference_fasta}")
SNP_file = config["SNP"]
sample_df = pd.read_csv('sample.tab', sep='\t', index_col=['SAMPLE'])
dict_Tn5_fastq   = {}
for sample_idx in sample_df.index:
	dict_Tn5_fastq.setdefault(sample_idx,{})
	with open(sample_df.at[sample_idx,"FASTQ_R1_FOFN"],'r') as fp:
		fq1_list = [x for x in fp.read().split('\n') if x.strip() != '']
	with open(sample_df.at[sample_idx,"FASTQ_R2_FOFN"],'r') as fp:
		fq2_list = [x for x in fp.read().split('\n') if x.strip() != '']
	if len(fq1_list) != len(fq2_list):
		print ("ERROR: # of files in fq1 and fq2 not matched")
		sys.exit()
	for i in range(len(fq1_list)):
		tmp_Tn5 = fq1_list[i].split('/')[-1].split('_')[0].split('.')[-1] #"split"+str(i+1) #Altos.GATCTATC_S1_R1_001.fastq.gz
		dict_Tn5_fastq[sample_idx].setdefault(tmp_Tn5,{}) # = [fq1_list[i],fq2_list[i]]
		dict_Tn5_fastq[sample_idx][tmp_Tn5].setdefault("R1",[])
		dict_Tn5_fastq[sample_idx][tmp_Tn5].setdefault("R2",[])
		dict_Tn5_fastq[sample_idx][tmp_Tn5]["R1"].append(fq1_list[i])
		dict_Tn5_fastq[sample_idx][tmp_Tn5]["R2"].append(fq2_list[i])

def find_ref(wildcards):
	return '/'.join(sample_df.at[wildcards.sample, "REF"].split('/')[0:-1])+'/'
def find_ref_idx(wildcards):
	return sample_df.at[wildcards.sample, "REF"] + '.fai'
def find_gtf_idx(wildcards):
	return sample_df.at[wildcards.sample, "GTF"]

def Tn5_fastq_R1(wildcards):
	return dict_Tn5_fastq[wildcards.sample][wildcards.Tn5]["R1"]
def Tn5_fastq_R2(wildcards):
	return dict_Tn5_fastq[wildcards.sample][wildcards.Tn5]["R2"]

rule all:
	input:
		expand("{sample}.merged.pass_cell_id.txt",sample=sample_df.index.get_level_values('SAMPLE')),
		##expand("scbs_data/scbs_filtered_{sample}.VMR.profile.csv",sample=sample_df.index.get_level_values('SAMPLE')),
		expand("allcools_file_list.{sample}.txt",sample=sample_df.index.get_level_values('SAMPLE')),
		#expand("tmp_bismark2extract_Tn5_path.{sample}.txt",sample=sample_df.index.get_level_values('SAMPLE')),

rule merge_Tn5_fastqs:
	input:
		R1_list = Tn5_fastq_R1,
		R2_list = Tn5_fastq_R2,
	output:
		combined_R1 = "Tn5_fastq/{sample}/{sample}_{Tn5}.R1.fq.gz",
		combined_R2 = "Tn5_fastq/{sample}/{sample}_{Tn5}.R2.fq.gz",
		Tn5_fastq_flag = "Tn5_fastq/{sample}/{sample}_{Tn5}.complete"
	threads: 1
	shell:
		'''
		cat {input.R1_list} > {output.combined_R1}
		cat {input.R2_list} > {output.combined_R2}
		touch Tn5_fastq/{wildcards.sample}/{wildcards.sample}_{wildcards.Tn5}.complete
		'''

rule trim_reads:
	input:
		R1 = rules.merge_Tn5_fastqs.output.combined_R1,
		R2 = rules.merge_Tn5_fastqs.output.combined_R2,
	output:
		trim_R1 = "trimmed_fastq/{sample}/{sample}_{Tn5}.trimmed.paired.R1.fq.gz",
		trim_R2 = "trimmed_fastq/{sample}/{sample}_{Tn5}.trimmed.paired.R2.fq.gz",
		tmp_fastq_flag = "tmp_fastq_dir/{sample}_{Tn5}/{sample}_{Tn5}.complete"
	threads: 4
	params:
		prefix_R1 = lambda wildcards, input: input.R1.split('/')[-1],
		prefix_R2 = lambda wildcards, input: input.R2.split('/')[-1]
	shell:
		'''
		mkdir -p tmp_fastq_dir/{wildcards.sample}_{wildcards.Tn5}
		cd tmp_fastq_dir/{wildcards.sample}_{wildcards.Tn5}
		perl {SDIR}/sciMETv2_workflows/sciMET_trim.pl -t {threads} -1 ../../{input.R1} -2 ../../{input.R2} -O ../../trimmed_fastq/{wildcards.sample}/{wildcards.sample}_{wildcards.Tn5}
		touch {wildcards.sample}_{wildcards.Tn5}.complete
		'''

rule bismark_align:
	input:
		BISMARK_REF = reference_dir,
		R1 = rules.trim_reads.output.trim_R1,
		R2 = rules.trim_reads.output.trim_R2,
	output:
		Tn5_nsrt = "tmp/{sample}/{sample}_{Tn5}.nsrt.bam",
		Tn5_align_report = "align_report_bismark/{sample}_{Tn5}.align_report.txt"
	threads: 8
	params:
		new_threads = 4,

	shell:
		'''
		perl {SDIR}/sciMETv2_workflows/sciMET_align_pe_threads.pl -t {params.new_threads} -R {input.BISMARK_REF} -1 {input.R1} -2 {input.R2} -O {wildcards.sample}_{wildcards.Tn5}
		mv {wildcards.sample}_{wildcards.Tn5}.nsrt.bam tmp/{wildcards.sample}/
		mv {wildcards.sample}_{wildcards.Tn5}.align_report.txt align_report_bismark/
		rm -rf {wildcards.sample}_{wildcards.Tn5}.*
		rm -rf {wildcards.sample}_{wildcards.Tn5}
		'''

rule remove_duplicate:
	input:
		Tn5_nsrt = "tmp/{sample}/{sample}_{Tn5}.nsrt.bam"
	output:
		compl = "temp_complexity/{sample}.{Tn5}.complexity.txt",
		rmdup_bam = "temp_rmdup_bam/{sample}.{Tn5}.bbrd.q10.bam",
	threads: 4
	shell:
		'''
		perl {SDIR}/sciMETv2_workflows/sciMET_rmdup_pe.pl -t {threads} -q {min_map_score} -O {wildcards.sample}.{wildcards.Tn5} {input.Tn5_nsrt}
		mv {wildcards.sample}.{wildcards.Tn5}.complexity.txt temp_complexity/
		mv {wildcards.sample}.{wildcards.Tn5}.bbrd.q10.bam temp_rmdup_bam/
		'''

rule namesort_Tn5_bam:
	input:
		Tn5_rmdup_bam = "temp_rmdup_bam/{sample}.{Tn5}.bbrd.q10.bam"
	output:
		Tn5_nsrt_rmdup_bam = temp("tmp2/{sample}/{sample}.{Tn5}.bbrd.q10.nsrt.bam"),
	threads: 4
	shell:
		'''
		samtools sort -n -@{threads} -o {output.Tn5_nsrt_rmdup_bam} {input.Tn5_rmdup_bam} 
		'''

rule splitbam_by_bc:
	input:
		bam = "tmp2/{sample}/{sample}.{Tn5}.bbrd.q10.nsrt.bam",
		complexity = "temp_complexity/{sample}.{Tn5}.complexity.txt"
	output:
		bc_split_complete = "filtered_bc_bam/{sample}_{Tn5}/split.complete",
	threads: 2
	shell:
		'''
		python {SDIR}/scripts/splitbam_bc.py {input.bam} {input.complexity} {output.bc_split_complete} {wildcards.sample} {num_uniq_reads}
		'''

rule rmdup_by_bc:
	input:
		bc_split_complete = "filtered_bc_bam/{sample}_{Tn5}/split.complete",
	output:
		bc_rmdup_complete = "filtered_bc_rmdup_bam/{sample}_{Tn5}/rmdup.complete",
	threads: 1
	shell:
		'''
		python {SDIR}/scripts/rmdup_bc.py {input.bc_split_complete} {output.bc_rmdup_complete}
		'''

rule bismark_extract_by_bc:
	input:
		bc_rmdup_complete = "filtered_bc_rmdup_bam/{sample}_{Tn5}/rmdup.complete",
	output:
		bismark_extract_bc_complete = "filtered_bc_bismark2extract/{sample}_{Tn5}/bismark_extract.complete",
	threads: 4
	shell:
		'''
		python {SDIR}/scripts/bismark_extract_bc.py {input.bc_rmdup_complete} {output.bismark_extract_bc_complete}
		'''

rule bismark_cytosine_report:
	input:
		bismark_extract_bc_complete = "filtered_bc_bismark2extract/{sample}_{Tn5}/bismark_extract.complete",
		ref_file = reference_fasta,
	output:
		bc_CXreport_complete = "filtered_bc_CXreport/{sample}_{Tn5}/coverage2cytosine.complete",
	threads: 2
	resources:
		load = 50
	shell:
		'''
		python {SDIR}/scripts/bismark_CXreport_bc.py {input.ref_file} {input.bismark_extract_bc_complete} {output.bc_CXreport_complete}
		'''

rule CXreport2allc:
	input:
		bc_CXreport_complete = "filtered_bc_CXreport/{sample}_{Tn5}/coverage2cytosine.complete",
		cytosine_SNP_file = SNP_file,
	output:
		CXreport2allc_complete = "filtered_bc_allcools/{sample}_{Tn5}/CXreport2allc.complete",
	threads: 1
	resources:
		load = 50
	shell:
		'''
		python {SDIR}/scripts/make_allc_file.py {input.cytosine_SNP_file} {input.bc_CXreport_complete} {output.CXreport2allc_complete} 
		'''

#rule sort_bc_bam:
#	input:
#		bc_split_complete = "filtered_bc_bam/{sample}_{Tn5}/split.complete",
#	output:
#		sort_bc_bam_complete = "filtered_bc_bam/sorted/{sample}_{Tn5}/sort.complete",
#	threads: 2
#	shell:
#		'''
#		python {SDIR}/scripts/sort_bc_bam.py {input.bc_split_complete} {output.sort_bc_bam_complete} 
#		'''
#
#rule allc_step1:
#	input:
#		sorted_bc_bam_Tn5_dir_flag = rules.sort_bc_bam.output.sort_bc_bam_complete,
#		ref_file = reference_fasta,
#	output:
#		bc_allc_step1_complete = "filtered_bc_allc/{sample}_{Tn5}/allc_step1.complete",
#	threads: 4
#	shell:
#		'''
#		python {SDIR}/scripts/allc_step1.py {input.ref_file} {input.sorted_bc_bam_Tn5_dir_flag} {output.bc_allc_step1_complete} {min_map_score} 
#		'''

#rule allc_step2_autosome_rmSNP:
#	input:
#		bc_allc_Tn5_dir_flag = rules.allc_step1.output.bc_allc_step1_complete,
#		cytosine_SNP_file = SNP_file,
#	output:
#		allc_step2_complete = "filtered_bc_allc/autosome_rmSNP/{sample}_{Tn5}/allc_step2.complete",
#	threads: 1
#	resources:
#		load = 50
#	shell:
#		'''
#		python {SDIR}/scripts/allc_step2_rmSNP.py {input.cytosine_SNP_file} {input.bc_allc_Tn5_dir_flag} {output.allc_step2_complete} 
#		'''

rule make_bismark2extract_list:
	input:
		bismark2extract_complete_flag = lambda wildcards: expand('filtered_bc_bismark2extract/{{sample}}_{Tn5}/bismark_extract.complete',Tn5=dict_Tn5_fastq[wildcards.sample].keys())
	output:
		bismark2extract_Tn5_path = "tmp_bismark2extract_Tn5_path.{sample}.txt",
	threads: 1
	shell:
		'''
		ls {input.bismark2extract_complete_flag} >| {output.bismark2extract_Tn5_path}
		'''

rule make_allc_table:
	input:
		allc_step2_complete_flag = lambda wildcards: expand('filtered_bc_allcools/{{sample}}_{Tn5}/CXreport2allc.complete',Tn5=dict_Tn5_fastq[wildcards.sample].keys())
	output:
		allc_Tn5_path = "tmp_allcools_Tn5_path.{sample}.txt",
		allc_file_list = "allcools_file_list.{sample}.txt",
	threads: 1
	shell:
		'''
		ls {input.allc_step2_complete_flag} >| {output.allc_Tn5_path}
		python {SDIR}/scripts/allc_file_path.py {output.allc_Tn5_path} {output.allc_file_list}
		'''

#################################################3
rule scbs_prepare:
	input:
		bismark_extract_bc_complete = lambda wildcards: expand('filtered_bc_bismark2extract/{{sample}}.{Tn5}/bismark_extract.complete',Tn5=dict_Tn5_fastq[wildcards.sample].keys())
	output:
		Tn5_bismark_extract_dir = "Tn5_bismark_extract_dir.{sample}.txt",
		scbs_prepare_out = "scbs_data/scbs_{sample}/run_info.txt",
	threads: 32
	shell:
		'''
		ls {input.bismark_extract_bc_complete} > {output.Tn5_bismark_extract_dir}
		python {SDIR}/scripts/scbs_prepare_run.py {output.Tn5_bismark_extract_dir} {output.scbs_prepare_out}
		'''

rule scbs_filter:
	input:
		rules.scbs_prepare.output.scbs_prepare_out,
	output:
		scbs_filter_out = "scbs_data/scbs_filtered_{sample}/run_info.txt",
	threads: 12
	shell:
		'''
		touch {input[0]}
		scbs filter --min-sites 200000 --max-sites 10000000 --min-meth 50 --max-meth 85 scbs_data/scbs_{wildcards.sample} scbs_data/scbs_filtered_{wildcards.sample}
		'''

rule scbs_smooth:
	input:
		rules.scbs_filter.output.scbs_filter_out,
	output:
		scbs_smooth_out = "scbs_data/scbs_filtered_{sample}/smoothed/check.complete",
	threads: 12
	params:
		outdir = "scbs_data/scbs_filtered_{sample}/smoothed"
	shell:
		'''
		touch {input[0]}
		scbs smooth scbs_data/scbs_filtered_{wildcards.sample}
		touch {params.outdir}/smoothed.complete
		'''

rule scbs_scan:
	input:
		rules.scbs_smooth.output.scbs_smooth_out,
	output:
		scbs_scan_out = "scbs_data/scbs_filtered_{sample}.VMR.bed",
	threads: 48
	params:
		var_per = 0.05
	shell:
		'''
		touch {input[0]}
		scbs scan --var-threshold {params.var_per} --threads {threads} scbs_data/scbs_filtered_{wildcards.sample} {output.scbs_scan_out}
		'''

rule scbs_profile:
	input:
		smooth = rules.scbs_smooth.output.scbs_smooth_out,
		VMR = rules.scbs_scan.output.scbs_scan_out,
	output:
		scbs_profile_out = "scbs_data/scbs_filtered_{sample}.VMR.profile.csv",
	threads: 12
	shell:
		'''
		touch {input.smooth}
		scbs profile {input.VMR} scbs_data/scbs_filtered_{wildcards.sample} {output.scbs_profile_out}
		'''

rule filter_cell_id:
	input:
		rules.remove_duplicate.output.compl,
	output:
		"pass_cell_id/{sample}.{Tn5}.pass_cell_id.txt",
	threads: 1
	shell:
		'''
		python {SDIR}/scripts/pass_filter_cellID.py {input} {output} {num_uniq_reads}
		'''

rule merge_pass_cell_id:
	input:
		pass_cell_file_list = lambda wildcards: expand('pass_cell_id/{{sample}}.{Tn5}.pass_cell_id.txt',Tn5=dict_Tn5_fastq[wildcards.sample].keys())
	output:
		merged_bam = "{sample}.merged.pass_cell_id.txt",
	threads: 1
	shell:
		'''
		cat {input.pass_cell_file_list} > {output.merged_bam}
		'''

rule methyl_extract:
	input:
		rules.remove_duplicate.output.rmdup_bam,
	output:
		CG_val= "{sample}.mCG.vals",
		CH_val= "{sample}.mCH.vals",
		CH_chrom_folder = "{sample}.CH.chrom_folders.txt",
		CG_chrom_folder = "{sample}.CG.chrom_folders.txt",
	threads: 16
	shell:
		'''
		perl {SDIR}/sciMETv2_workflows/sciMET_extract.pl -t {threads} -X -O {wildcards.sample} {input}
		'''

rule sort_bed:
	input:
		CG_chrom_folder = rules.methyl_extract.output.CG_chrom_folder,
		CH_chrom_folder = rules.methyl_extract.output.CH_chrom_folder,
	output:
		CG_chrom_sort_flag = "{sample}.CG.chrom_folders.sort.done",
		CH_chrom_sort_flag = "{sample}.CH.chrom_folders.sort.done",
	threads: 32
	shell:
		'''
		python {SDIR}/scripts/sort_chroms.py {input.CG_chrom_folder} {threads} {output.CG_chrom_sort_flag}
		python {SDIR}/scripts/sort_chroms.py {input.CH_chrom_folder} {threads} {output.CH_chrom_sort_flag}
		'''

rule make_windows:
	input:
		fasta_idx = reference_fasta_fai,#find_ref_idx,
		mock_input = "{sample}.pass_cell_id.txt"
	output:
		window_file = "{sample}.windows.bed"
	threads: 4
	params:
		window_size = 50000,
		min_chrom = 10000000,
	shell:
		'''
		echo {input.mock_input}
		cat {input.fasta_idx} | cut -f 1,2 | awk "{{OFS="\\t"}};{{if (\$2 > {params.min_chrom}) print \$0}}" | bedtools makewindows -w {params.window_size} -g - > {output.window_file}
		'''

rule methyl_to_matrix:
	input:
		chrom_folder = "{sample}.{context}.chrom_folders.txt",
		win_file = "{sample}.windows.bed",
		pass_cell_file = "{sample}.pass_cell_id.txt",
	output:
		mtx_ratio = "{sample}.{context}.ratio.mtx",
		mtx_cov = "{sample}.{context}.cov.mtx",
		mtx_mat = "{sample}.{context}.mtx",
	threads: 12
	shell: # will change the code to process for each chrom
		'''
		echo {output.mtx_ratio}
		perl {SDIR}/sciMETv2_workflows/sciMET_meth2mtx.pl -F {input.chrom_folder} -B {input.win_file} -O {wildcards.sample} -C {input.pass_cell_file}
		'''

#rule make_sliding_windows:
#	input:
#		gtf = gtf_file,
#	output:
#		protein_coding_bed = "protein_coding.bed",
#		sliding_window_file = "protein_coding.sliding_windows.bed",
#	threads: 4
#	params:
#		max_feature_len = 100000,
#		upstream_len = 5000,
#		upstream_win_num = 50,
#		downstream_len = 5000,
#		downstream_win_num = 50,
#	shell:
#		'''
#		cat {input.gtf} | awk -v OFS="\\t" "{{if (\$3 == "gene") print \$0}}" | grep "protein_coding" | awk -F"\\t" "BEGIN {{OFS = FS}};{{print \$1,\$4-1,\$5,\$9,\$7}}" | grep "chr" >| {output.protein_coding_bed}
#		perl {SDIR}/sciMETv2_workflows/sciMET_featuresToScanBed.pl -M {params.max_feature_len} -S {params.upstream_len} -s {params.upstream_win_num} -E {params.downstream_len} -e {params.downstream_win_num} {output.protein_coding_bed} {output.sliding_window_file}
#		'''

#rule win_Meth:
#	input:
#		chrom_folder = "{sample}.{context}.chrom_folders.txt",
#		slidingwin_file = sliding_win_file,#"protein_coding.sliding_windows.bed",
#	output:
#		win_Meth_out = "{sample}_{context}_winMeth.{context}.txt",
#	threads: 4
#	shell: # will change the code to process for each chrom
#		'''
#		perl {SDIR}/sciMETv2_workflows/sciMET_getWindowMeth.pl -F {input.chrom_folder} -O {wildcards.sample}_{wildcards.context}_winMeth -B {input.slidingwin_file}
#		'''
#perl ../scimetv2_pipeline/sciMETv2/sciMET_getWindowMeth.pl -F kidney-A.CG.chrom_folders.txt -O kidney-A_CG_winMeth -B protein_coding.slidingwin.bed	
#perl ../scimetv2_pipeline/sciMETv2/sciMET_featuresToScanBed.pl -M 100000 -S 5000 -s 50 -E 5000 -e 50 ../ref/annotation/protein_coding.bed protein_coding.slidingwin.bed
