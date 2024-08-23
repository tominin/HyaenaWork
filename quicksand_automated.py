import sys
import pandas as pd
import csv
import re
import numpy as np
from scipy.stats import binom
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
import statsmodels.stats.proportion as smp
def main():
	pd.options.mode.chained_assignment = None
	input_all='20240618_data_for_tom_DC_MC_SE-Profile_AA75_AllSamples_quicksandFilter.csv'
	xdata_all=pd.read_csv(input_all)    #20240618_data_for_tom_DC_MC_SE-Profile_AA75_AllSamples_quicksandFilter.csv
	input_hom='hom_curated_input19062020.txt'
	xdata6=pd.read_table(input_hom)       #hom_curated_input19062020.txt
	layers='layers_ZJ_04012021.txt'
	layers_read=pd.read_table(layers)       #layers_ZJ_04012021.txt
	what_chamber=sys.argv[1]
	screen_info=[]
	what_screen=xdata6['Screen'].tolist()
	for item in what_screen:
		if item == 'First':
			screen_info.append(True)
		else:
			screen_info.append(False)
	mgmldata_all_cols=xdata6[screen_info]
	keepers=[0, 4, 11]    #SampleID, IndexLibID, chamber
	mgmldata=mgmldata_all_cols.iloc[:, keepers]
	layerdata=pd.merge(layers_read, mgmldata, on='SampleID')   #SampleID, marker, layer, IndexLibID, chamber
	layerdata['Layer'] = layerdata['Layer'].astype(str)
	layer_list=layerdata['Layer'].tolist()
	new_layer_list=[]
	for item in layer_list:
		fh_g=item.split('.')
		if re.search('[\+|\/|dMP|pdd]', fh_g[0]):
			new_layer_list.append(0)
		else:
			holder=fh_g[0].strip(' ?')
			holder=str(holder)
			new_layer_list.append(holder)
	layerdata = layerdata.drop('Layer', axis='columns')
	layerdata['Layer'] = new_layer_list
	input3='layerdates_ZJ.txt'
	xdata3 = pd.read_table(input3)   #layerdates_ZJ2.txt
	xdata3 = xdata3[xdata3['Chamber'] == what_chamber]   #only MAIN
	file_path1='hgA_count_postfilt.txt'
	file_path2='hgB_count_postfilt.txt'
	file_path3='hgD_count_postfilt.txt'
	file_path4='hgC_count_postfilt.txt'
	col_names=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
	def skip_comments(file):
		splicer = []
		with open(file, 'r') as fh_temp:
			data_lines=fh_temp.readlines()
			for line in data_lines:
				if re.search('#', line):
					splicer.append(False)
				else:
					splicer.append(True)
			data_struc=pd.read_table(file, header=None, names=col_names)
			data_fin=data_struc[splicer]
		return (data_fin)
	haploA=skip_comments(file_path1)
	haploB=skip_comments(file_path2)
	haploD=skip_comments(file_path3)
	haploC=skip_comments(file_path4)
	def binom_conf_int(successes, total, conf_level=0.95):
		return smp.proportion_confint(successes, total, alpha=1-conf_level, method='binom_test')
	def process_dps(dpfile, bamend, family, group, stats):
		xdata_fam=xdata_all[xdata_all['Family'] == family]   #Only keeps Hyaenidae family
		dptest=dpfile
		dptest_edit=dptest.replace(to_replace=bamend, value = "", regex=True)    #removes bam file suffix leaving only CapLibIDCoreDB marker
		dptest_fin=dptest_edit.rename(columns={'a': 'CapLibIDCoreDB'})
		haplo=pd.merge(xdata_fam, dptest_fin, on='CapLibIDCoreDB')
		col_index1=54    #opp hits
		col_index2=55    #good hits
		tempa=haplo.iloc[:, col_index1].tolist()
		tempb=haplo.iloc[:, col_index2].tolist()
		splicer=[]
		for j, k in zip(tempa, tempb):
			j=int(j)
			k=int(k)
			if k<3:    #removes any samples that had lower than 3 hits
				k=0
				j=0
			total_runs=j+k
			if total_runs>0:
				lower_bound, upper_bound=stats(k, total_runs)   #binomial test
				x=lower_bound*100
				if x >= 10:
					splicer.append(True)
				else:
					splicer.append(False)
			else:
				splicer.append(False)
		test=haplo[splicer]
		keeperz=[0, 1, 8, 51, 58]
		sigfam=test.iloc[:, keeperz]   #SampleID, IndexLibID, CapLibID, perc(of good hits), marker/profile
		if not sigfam.empty:
			sigfam.loc[:, 'group']=group   #adds a new column with the specific haplogroup
			sigfam.rename(columns={'f' : 'perc'}, inplace=True)
		finalsig=pd.merge(sigfam, layerdata, on='SampleID')
		return(finalsig)
	haploA_fin=process_dps(haploA, '(\\.Hyaenidae.Crocuta_crocuta_deduped_bedfiltered.bam)', 'Hyaenidae', 'HaploA', binom_conf_int)
	haploB_fin=process_dps(haploB, '(\\.Hyaenidae.Crocuta_crocuta_deduped_bedfiltered.bam)', 'Hyaenidae', 'HaploB', binom_conf_int)
	haploD_fin=process_dps(haploD, '(\\.Hyaenidae.Crocuta_crocuta_deduped_bedfiltered.bam)', 'Hyaenidae', 'HaploD', binom_conf_int)
	haploC_fin=process_dps(haploC, '(\\.Hyaenidae.Crocuta_crocuta_deduped_bedfiltered.bam)', 'Hyaenidae', 'HaploC', binom_conf_int)
	haphyae=pd.concat([haploA_fin, haploB_fin, haploD_fin, haploC_fin], ignore_index=True)
	def process_sup_layer_ch(gen_mam):
		results=[]
		gengroup = sorted(gen_mam['group'].unique())  #sorted list of haplogroup names
		genlayer = gen_mam['Layer'].unique()
		genlayer = [int(i) for i in genlayer]
		genlayer.sort()
		genlayer = [str(i) for i in genlayer]    #sorted list of each layers
		gen_mam = gen_mam.drop_duplicates(subset=['SampleID', 'group'])   #subset powerful tool!! would save me alot of code in future
		samplesperlayer = gen_mam['Layer'].value_counts().reset_index()   #counts all the samples for each layer from all groups
		samplesperlayer.columns = ['Layer', 'Freq']
		for sfam in gengroup:
			for s1 in genlayer:
				sum_set = gen_mam[(gen_mam['Layer'] == s1) & (gen_mam['group'] == sfam)]
				samps_row = samplesperlayer[samplesperlayer['Layer'] == s1]
				if not samps_row.empty:
					samps = samps_row['Freq'].values[0]
				else:
					samps = 0    #samps = total count of samples for all 3 haplogroups
				total_group = sum_set['SampleID'].nunique()   #number of uniq sampleid's for haplo group in a specific layer
				por_perlayer = (total_group / samps * 100) if samps != 0 else 0
				results.append({
			'Chamber': what_chamber,
			'layer': s1,
			'sfam': sfam,
			'por_perlayer': por_perlayer,
			'samps': samps
						})
		g = pd.DataFrame(results)
		layerg = pd.merge (g, xdata3, left_on=['layer'], right_on=['layer'])
		layerg['center'] = (layerg['start'] + layerg['end']) / 2
		return(layerg)
	porp_hyae = process_sup_layer_ch(haphyae)
	fig, ax = plt.subplots(figsize=(12, 8))
	indent={}
	colors = {'HaploA': '#8b4513', 'HaploB': '#53290b', 'HaploD': '#c5a289', 'HaploC': '#FE420F'}
	for _, row in porp_hyae[porp_hyae['sfam'] == 'HaploB'].iterrows():
		ax.broken_barh([(0, 100)], (row['start'], row['end'] - row['start']), facecolors=colors['HaploB'], label='barh')
	for _, row in porp_hyae[porp_hyae['sfam'] == 'HaploA'].iterrows():
		ax.broken_barh([(100 - row['por_perlayer'], row['por_perlayer'])], (row['start'], row['end'] - row['start']), facecolors=colors['HaploA'], label='barh')
	for _, row in porp_hyae[porp_hyae['sfam'] == 'HaploD'].iterrows():
		ax.broken_barh([(0, row['por_perlayer'])], (row['start'], row['end'] - row['start']), facecolors=colors['HaploD'], label='barh')
		indent[row['start']] = row['por_perlayer']
	for _, row in porp_hyae[porp_hyae['sfam'] == 'HaploC'].iterrows():
		for key, value in indent.items():
			if row['start'] == key:
				ax.broken_barh([(value, (row['por_perlayer']))], (row['start'], row['end'] - row['start']), facecolors=colors['HaploC'], label= 'barh')
	ax.scatter(porp_hyae['samps'], porp_hyae['center'], color='white', edgecolor='black', zorder=5)
	ax.set_xlabel("Percent")
	ax.set_ylabel("Time (ka)")
	ax.set_title("Cave Hyaena Chronological Distrubution")
	ax.invert_yaxis()
	ax.set_ylim(300, 0)
	ax.set_yticks(range(0, 301, 50))
	secax = ax.secondary_xaxis('top')
	secax.set_xlabel('Total no. of samples')
	box=ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	legend_elements=[Patch(facecolor=colors['HaploA'], label='Haplogroup A'),
			Patch(facecolor=colors['HaploB'], label='Haplogroup B'),
			Patch(facecolor=colors['HaploC'], label='Haplogroup C'),
			Patch(facecolor=colors['HaploD'], label='Haplogroup D')]
	ax.legend(handles=legend_elements, bbox_to_anchor=(1, 1), loc='upper left')
	plt.savefig("denisova_maincham_hyae.png")
main()

