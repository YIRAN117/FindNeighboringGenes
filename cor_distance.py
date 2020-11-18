#goal: check whether transcription correlation is affected by gene distance
#function: generate cor and distance curve lines for gene pairs to check 
#input: correlation table between gene pairs; gene gff
#output: distance cor table file and figure

import re
import os

#read in gene pair information including pearson correlation, p-value, and gff
os.system('''head -n 4335 /Users/yiranli/Dropbox/lncRNA_article/data/r >/Users/yiranli/Dropbox/lncRNA_article/data/r1''')
g_input_gff_file="/Users/yiranli/Dropbox/lncRNA_article/lncRNA_checked_backed_merged.gff"
g_gene_info_dict={}
g_input_cor_file="/Users/yiranli/Dropbox/lncRNA_article/data/r1"
g_input_pvalue_file="/Users/yiranli/Dropbox/lncRNA_article/data/fdr"
g_output_lncRNA_mRNA="/Users/yiranli/Dropbox/lncRNA_article/data/cor/lncRNA_mRNA_cor_distance"
g_output_mRNA_mRNA="/Users/yiranli/Dropbox/lncRNA_article/data/cor/mRNA_mRNA_cor_distance"

#build gene information dict
def gene_info(gff_file):
    with open(g_input_gff_file) as f1:
    	a_dict={}
    	for line in f1:
    		if not line.startswith("#"):
	            features=line.rstrip().split("\t")
	            if features[2]=="exon" or features[2]=="mRNA":
	                chro=features[0]
	                gene_id=features[8].split(";")[0].replace("ID=","").replace("Cgd","cgd")
	                if gene_id.startswith("Parent="):
	                	gene_id=re.search(r'Parent=.*',gene_id).group().replace("Parent=","")+"-RA"
	                strand=features[6]
	                if strand=="+":
	                	start=int(features[3]);end=int(features[4])
	                else:
	                	start=int(features[4]);end=int(features[3])
	                if gene_id not in g_gene_info_dict.keys():
	                	a_dict[gene_id]={"chro":chro,"strand":strand,"start":start,"end":end}
	                else:
	                	pass
    return a_dict

#calculate gene pair distance
def gene_distance(gene1_dict,gene2_dict):
	if gene1_dict["chro"]==gene2_dict["chro"]:
		if len(set(range(min(gene1_dict["start"],gene1_dict["end"]),max(gene1_dict["start"],gene1_dict["end"])+1)).intersection(set(range(min(gene2_dict["start"],gene2_dict["end"]),max(gene2_dict["start"],gene2_dict["end"])+1))))>2:
			d=0#overlap
		else:
			d=min(abs(gene1_dict["start"]-gene2_dict["start"]),abs(gene1_dict["start"]-gene2_dict["end"]),abs(gene1_dict["end"]-gene2_dict["end"]),abs(gene1_dict["end"]-gene2_dict["start"]))
	else:
		d=-1
	return d

#add gene expression correlation information
def cor_info(a_cor_table):
	gene_order_dict={}
	a_dict={}
	with open(g_input_cor_file) as f2:
		header=f2.readline().replace("Cgd","cgd")
		gene_order=header.rstrip().split("\t")
		gene_order_len=len(gene_order)
		print(gene_order_len)
		#gene_order=[x.replace("mRNA_","") for x in gene_order]
		gene_order_itor=enumerate(gene_order)#[(0,gene1),(1,gene2)]
		for i in gene_order_itor:
			gene_order_dict[i[1]]=i[0] #{mRNA_Cgd2_140-RA:0,mRNA_Cgd2_1810-RA :1,mRNA_Cgd2_2990-RA:2}
		for line in f2:
			features=line.rstrip().split("\t")
			gene_id=features[0].replace("Cgd","cgd")
			cor_list=features[1:]
			#print(len(cor_list))	
			for i in range(gene_order_dict[gene_id]+1,gene_order_len):
				key=tuple(sorted([gene_id,gene_order[i]]))
				try:
					a_dict[key]=float(cor_list[i])
				except:
					a_dict[key]=10
	return a_dict

#add p-value information
def pvalue_info(a_pvalue_table):
	with open(g_input_pvalue_file) as f3:
		gene_order_dict={}		
		a_dict={}
		header=f3.readline().replace("Cgd","cgd")
		gene_order=header.rstrip().split("\t")
		gene_order_len=len(gene_order)
		#gene_order=[x.replace("mRNA_","") for x in gene_order]
		gene_order_itor=enumerate(gene_order)
		for i in gene_order_itor:
			gene_order_dict[i[1]]=i[0] #{mRNA_Cgd2_140-RA:0,mRNA_Cgd2_1810-RA :1,mRNA_Cgd2_2990-RA:2}
		for line in f3:
			features=line.rstrip().split("\t")
			gene_id=features[0].replace("Cgd","cgd")
			pvalue_list=features[1:]	
			for i in range(gene_order_dict[gene_id]+1,gene_order_len):
				key=tuple(sorted([gene_id,gene_order[i]]))
				try:
					a_dict[key]=float(pvalue_list[i])
				except:
					a_dict[key]=10			
	return a_dict


g_gene_info_dict=gene_info(g_input_gff_file)
g_cor_dict=cor_info(g_input_cor_file)
g_pvalue_dict=pvalue_info(g_input_pvalue_file)

with open(g_output_mRNA_mRNA,"w") as write1:
	write1.write('gene\tdistance\tcorrelation\tfdr\n')
	with open(g_output_lncRNA_mRNA,"w") as write2:
		write2.write('gene\tdistance\tcorrelation\tfdr\n')
		for key,value in g_cor_dict.items():
			if g_pvalue_dict[key] < 2:# and abs(g_cor_dict[key]) >0.5:
				if "mRNA" in key[1] and "mRNA" not in key[0]:
					try:
						write2.write("{}:{}\t{}\t{}\t{}\n".format(key[0],key[1],gene_distance(g_gene_info_dict[key[0]],g_gene_info_dict[key[1].replace("mRNA_","")]),value,g_pvalue_dict[key]))
					except:
						pass
						print(key)
				if "mRNA" in key[1] and "mRNA" in key[0]:
					try:
						write1.write("{}:{}\t{}\t{}\t{}\n".format(key[0],key[1],gene_distance(g_gene_info_dict[key[0].replace("mRNA_","")],g_gene_info_dict[key[1].replace("mRNA_","")]),value,g_pvalue_dict[key]))
					except:
						pass
						print(key)

#make figure to check whether transcription correlation is affected by gene distance
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt

plt.style.use('classic')
import pandas as pd
lncRNA_mRNA_table=pd.read_csv(g_output_lncRNA_mRNA,sep="\t")
ax = sns.regplot(x="distance", y="correlation", data=lncRNA_mRNA_table)
plt.show()






                







