#!/usr/bin/python3
import sys
import matplotlib.pyplot as plt
import csv
import pandas as pd
import re
import numpy as np
import math
def main():
	for line in sys.stdin:
		in_path=line.rstrip()
		break
	dp_path='/mnt/scratch/thomas_harris_snell/hyena_mtDNA/hyena_diagnostic_positions/non_filtered_diag_pos/tom/'
	outlier_dict={}
	inlier_dict={}
	a_file=open(dp_path + '/dp_hgA.txt', 'r')   #DIAGNOSTIC POS FILE
	a_fh=pd.read_table(a_file, sep='\t')
	b_file=open(dp_path + '/dp_hgB.txt', 'r')
	b_fh=pd.read_table(b_file, sep='\t')
	c_file=open(dp_path + '/dp_hgC.txt', 'r')
	c_fh=pd.read_table(c_file, sep='\t')
	d_file=open(dp_path + '/dp_hgD.txt', 'r')
	d_fh=pd.read_table(d_file, sep='\t')
	a_concat_file=open(in_path + '/hgA_prefilt.txt', 'r')
	b_concat_file=open(in_path + '/hgB_prefilt.txt', 'r')
	c_concat_file=open(in_path + '/hgC_prefilt.txt', 'r')
	d_concat_file=open(in_path + '/hgD_prefilt.txt', 'r')
	a_splicer=[]
	b_splicer=[]
	c_splicer=[]
	d_splicer=[]
	def file1sorter(fh, splicer):
		pos_list=[]
		keeper=0
		temp=fh.iloc[:, keeper]
		pre_filt=temp.tolist()
		for item in pre_filt:
			if item[0] != '#':
				int_item=int(item)
				pos_list.append(int_item)
				splicer.append(True)
			else:
				splicer.append(False)
		return(pos_list, splicer)
	a_dp_data=file1sorter(a_fh, a_splicer)
	a_spliced=a_fh[a_splicer]
	b_dp_data=file1sorter(b_fh, b_splicer)
	b_spliced=b_fh[b_splicer]
	c_dp_data=file1sorter(c_fh, c_splicer)
	c_spliced=c_fh[c_splicer]
	d_dp_data=file1sorter(d_fh, d_splicer)
	d_spliced=d_fh[d_splicer]
	def file2sorter(file2):
		list=[]
		data2=file2.readlines()
		data3=data2[1:]
		for line in data3:
			f=line.split("\t")
			temp=f[7]
			temp2=f[8]
			temp3=f[9]
			no_symbol=re.split(r"[:|,]", temp)
			no_symbol2=re.split(r"[:|,]", temp2)
			no_symbol3=re.split(r"[:|,]", temp3)
			cat_no_symbol=no_symbol + no_symbol2 + no_symbol3
			for item in cat_no_symbol:
				if item.isdigit():
					temp2=int(item)
					list.append(temp2)
				else:
					continue
		return(list)
	a_count_data=file2sorter(a_concat_file)
	b_count_data=file2sorter(b_concat_file)
	c_count_data=file2sorter(c_concat_file)
	d_count_data=file2sorter(d_concat_file)
	def dictmaker(dp_data, count_data):
		pos_count={}
		for diag_pos in dp_data:
			coverage=count_data.count(diag_pos)
			pos_count[diag_pos]=coverage
		return(pos_count)
	a_dict=dictmaker(a_dp_data[0], a_count_data)
	b_dict=dictmaker(b_dp_data[0], b_count_data)
	c_dict=dictmaker(c_dp_data[0], c_count_data)
	d_dict=dictmaker(d_dp_data[0], d_count_data)
	def maths_time(ruler_dict):
		total=len(ruler_dict)
		cumulative=0
		sum_var=0
		count=0
		for key, item in ruler_dict.items():   #key=pos item=coverage
			cumulative=cumulative+item
		mean=cumulative/total
		for key2, item2 in ruler_dict.items():
			variation=item2-mean
			count=count+1
			var_square=variation*variation
			sum_var=sum_var+var_square
		penul=sum_var/count
		std_dev=math.sqrt(penul)
		threshold=mean+(2*std_dev)
		return(threshold)			#seperate positions that are > 2xthe mean value of coverage
	a_std_var=maths_time(a_dict)
	b_std_var=maths_time(b_dict)
	c_std_var=maths_time(c_dict)
	d_std_var=maths_time(d_dict)
	def remover_of_insignificant(ruler_dict, std_var):
		outlier_dict={}
		inlier_dict={}
		for key, item in ruler_dict.items():
			key=str(key)
			if item > std_var:
				outlier_dict[key]=item
			else:
				inlier_dict[key]=item
		out=list(outlier_dict.keys())
		out.sort()
		return (out)
	a_outliers=remover_of_insignificant(a_dict, a_std_var)
	a_outlierz=', '.join(a_outliers)
	print ('Haplo As outlying diag pos:' + a_outlierz)
	b_outliers=remover_of_insignificant(b_dict, b_std_var)
	b_outlierz=', '.join(b_outliers)
	print ('Haplo Bs outlying diag pos:' + b_outlierz)
	c_outliers=remover_of_insignificant(c_dict, c_std_var)
	c_outlierz=', '.join(c_outliers)
	print ('Haplo Cs outlying diag pos:' + c_outlierz)
	d_outliers=remover_of_insignificant(d_dict, d_std_var)
	d_outlierz=', '.join(d_outliers)
	print ('Haplo Ds outlying diag pos:' + d_outlierz)
	def df_output(spliced, out):
		splicer2=[]
		out2=[]
		keeper = 0
		temp = spliced.iloc[:, keeper]
		pre_filt = temp.tolist()
		for item in out:
			new=str(item)
			out2.append(new)
		fin_df = spliced[~spliced['##position'].isin (out2)]
		return(fin_df)
	a_df=df_output(a_spliced, a_outliers)
	b_df=df_output(b_spliced, b_outliers)
	c_df=df_output(c_spliced, c_outliers)
	d_df=df_output(d_spliced, d_outliers)
	out_filepath=sys.argv[1]  #mnt/scratch/thomas_harris_snell/hyena_mtDNA/hyena_diagno>
	a_df.to_csv(out_filepath + '/filt_dp_hgA.txt', index=None, sep='\t', mode='w')
	b_df.to_csv(out_filepath + '/filt_dp_hgB.txt', index=None, sep='\t', mode='w')
	c_df.to_csv(out_filepath + '/filt_dp_hgC.txt', index=None, sep='\t', mode='w')
	d_df.to_csv(out_filepath + '/filt_dp_hgD.txt', index=None, sep='\t', mode='w')
	text_file = open('del_me2.txt', 'w')
	text_file.write(out_filepath)
	text_file.close()
main()


