# function:find lncRNA homology between two species based on relative position to mRNAs. 
# input:[lncrna gff, mRNA gff,mRNA OG file]-for each species, ortholog file,  
# output: conserved cpar lncrna and related mrna

from pybedtools import BedTool
import re
#global variable

pf3d7_gene_OG_file="/Users/yiranli/Dropbox/master_assembly/homology/pf3d7_gene_OG_info.tab"
cpar_gene_OG_file="/Users/yiranli/Dropbox/master_assembly/homology/cpar_gene_OG_info.tab"
pf3d7_lncRNA_info_file="/Users/yiranli/Dropbox/master_assembly/homology/2015_pf3d7_lncRNA.tab"
pf3d7_gene_info_file="/Users/yiranli/Dropbox/master_assembly/homology/PlasmoDB-45_Pfalciparum3D7.gff"
cpar_gene_infor_file="/Users/yiranli/Dropbox/master_assembly/homology/CryptoDB-45_CparvumIowaII.gff"
cpar_lncRNA_info_file="/Users/yiranli/Dropbox/master_assembly/list_clean/clean.gff"
cpar_lncRNA_list=[]
pf3d7_lncRNA_list=[]
#make lnc db
g_2015_pf3d7_lnc_db=""
g_clean_cpar_lnc_db=""
with open(pf3d7_lncRNA_info_file) as f1:
	header=f1.readline()
	entry_all=[]
	for line in f1:
		features=line.rstrip().split()
		lnc_id=features[0]
		pf3d7_lncRNA_list.append(lnc_id)
		chrom=features[1]
		start=features[2]
		end=features[3]
		strand=features[4]
		entry='{}\t{}\t{}\t{}\t1\t{}\n'.format(chrom,start,end,lnc_id,strand)
		entry_all.append(entry)
	bed="\n".join(entry_all)
	g_2015_pf3d7_lnc_db = BedTool(bed,from_string=True).sort()

with open(cpar_lncRNA_info_file) as f2:
	entry_all=[]
	for line in f2:
		features=line.rstrip().split("\t")
		try:
			lnc_id=re.search("TCONS_\d+",features[8]).group()
		except:
			lnc_id=re.search("MSTRG.\d+",features[8]).group()
		cpar_lncRNA_list.append(lnc_id)
		chrom=features[0]
		start=features[3]
		end=features[4]
		strand=features[6]
		entry='{}\t{}\t{}\t{}\t1\t{}\n'.format(chrom,start,end,lnc_id,strand)
		entry_all.append(entry)
	bed="\n".join(entry_all)
	g_clean_cpar_lnc_db = BedTool(bed,from_string=True).sort()


#make gene db
def make_gene_db(file):
	out_db={}
	with open(file) as handle:
		entry_all=[]
		for line in handle:
			features=line.rstrip().split("\t")
			if len(features)==9 and features[2]=="gene":
				chrom=features[0]
				start=features[3]
				end=features[4]
				strand=features[6]
				gene_id=re.search("ID=.*?;",features[8]).group().replace("ID=","").replace(";","")
				entry='{}\t{}\t{}\t{}\t1\t{}\n'.format(chrom,start,end,gene_id,strand)
				entry_all.append(entry)
		bed="\n".join(entry_all)
		out_db = BedTool(bed,from_string=True).sort()
	return out_db

g_2015_pf3d7_gene_db=make_gene_db(pf3d7_gene_info_file)
g_2015_cpar_gene_db=make_gene_db(cpar_gene_infor_file)

#make lnc and gene relationship 
g_pf3d7_upstream=g_2015_pf3d7_lnc_db.closest(g_2015_pf3d7_gene_db,io=True,id=True,D="a")
g_pf3d7_downstream=g_2015_pf3d7_lnc_db.closest(g_2015_pf3d7_gene_db,io=True,iu=True,D="a")
g_pf3d7_antisense=g_2015_pf3d7_lnc_db.intersect(g_2015_pf3d7_gene_db,wao=True,S=True)

g_cpar_upstream=g_clean_cpar_lnc_db.closest(g_2015_cpar_gene_db,io=True,id=True,D="a")
g_cpar_downstream=g_clean_cpar_lnc_db.closest(g_2015_cpar_gene_db,io=True,iu=True,D="a")
g_cpar_antisense=g_clean_cpar_lnc_db.intersect(g_2015_cpar_gene_db,wao=True,S=True)
#put relationship into dict
def make_relation_dict(input_db):
	out_dict={}
	for i in input_db:
		if i[-1]!=0 and i[-1]!=-1:
			out_dict[i[9]]=i[3]
	return out_dict 

def make_relation_dict_forward(input_db):
	out_dict={}
	for i in input_db:
		if i[-1]!=0 and i[-1]!=-1:
			out_dict[i[3]]=i[9]
	return out_dict 
g_pf3d7_dict_upstream={}
g_pf3d7_dict_downstream={}
g_pf3d7_dict_antisense={}
g_cpar_dict_upstream={}
g_cpar_dict_downstream={}
g_cpar_dict_antisense={}

g_pf3d7_dict_upstream=make_relation_dict_forward(g_pf3d7_upstream)
g_pf3d7_dict_downstream=make_relation_dict_forward(g_pf3d7_downstream)
g_pf3d7_dict_antisense=make_relation_dict(g_pf3d7_antisense)
g_cpar_dict_upstream=make_relation_dict_forward(g_cpar_upstream)
g_cpar_dict_downstream=make_relation_dict_forward(g_cpar_downstream)
g_cpar_dict_antisense=make_relation_dict_forward(g_cpar_antisense)


#make plasmo and cpar OG and mrna relation
g_pf3d7_OG_gene_dict={}
g_cpar_OG_gene_dict={}

def make_OG_gene_dict(file):
	out_dict={}
	with open(file) as handle:
		header=handle.readline()
		for line in handle:
			features=line.rstrip().split()
			gene=features[0]
			group=features[3]
			if group not in out_dict.keys():
				out_dict[group]=[gene]
			else:
				out_dict[group].append(gene)
	return out_dict

def make_OG_gene_dict_forward(file):
	out_dict={}
	with open(file) as handle:
		header=handle.readline()
		for line in handle:
			features=line.rstrip().split()
			gene=features[0]
			group=features[3]
			out_dict[gene]=group
	return out_dict
g_pf3d7_OG_gene_dict=make_OG_gene_dict(pf3d7_gene_OG_file)
g_pf3d7_OG_gene_dict_forward=make_OG_gene_dict_forward(pf3d7_gene_OG_file)
g_cpar_OG_gene_dict=make_OG_gene_dict_forward(cpar_gene_OG_file)

level2=[]
for i_lnc_cpar in cpar_lncRNA_list:
	try:#not exist because of cpar intergenic lncrna
		i_cpar_antisense=g_cpar_dict_antisense[i_lnc_cpar]
		i_cpar_antisense_OG=g_cpar_OG_gene_dict[i_cpar_antisense]
		i_pf3d7_antisense_OG_gene_list=g_pf3d7_OG_gene_dict[i_cpar_antisense_OG]
		for i in i_pf3d7_antisense_OG_gene_list:
			if i in g_pf3d7_dict_antisense.keys():
				level2.append(i_lnc_cpar)
				print("{}\t{}\t{}\tantisense ortholog exist in pf37d".format(i_lnc_cpar,i_cpar_antisense,i))
	except:
		pass

g_cpar_dict_neighbor={}
for i_lnc_cpar in cpar_lncRNA_list:
	try:
		upstream=g_cpar_dict_upstream[i_lnc_cpar]
		upstream_OG=g_cpar_OG_gene_dict[upstream]
		downstream=g_cpar_dict_downstream[i_lnc_cpar]
		downstream_OG=g_cpar_OG_gene_dict[downstream]
		value=[upstream_OG,downstream_OG]
		value.sort()
		g_cpar_dict_neighbor[i_lnc_cpar]=value
	except:
		pass


g_pf3d7_dict_neighbor={}
for i_lnc_pf3d7 in pf3d7_lncRNA_list:
	try:
		upstream=g_pf3d7_dict_upstream[i_lnc_pf3d7]
		upstream_OG=g_pf3d7_OG_gene_dict_forward[upstream]
		downstream=g_pf3d7_dict_downstream[i_lnc_pf3d7]
		downstream_OG=g_pf3d7_OG_gene_dict_forward[downstream]
		key=[upstream_OG,downstream_OG]
		key.sort()
		key=tuple(key)
		g_pf3d7_dict_neighbor[key]=i_lnc_pf3d7
	except:
		pass

for i in cpar_lncRNA_list:
	try:
		cpar_neighbor=g_cpar_dict_neighbor[i]
		cpar_neighbor=tuple(cpar_neighbor)
		if cpar_neighbor in g_pf3d7_dict_neighbor.keys():
			print('{} has level B ortholog:{}'.format(i,g_pf3d7_dict_neighbor[cpar_neighbor]))
	except:
		pass
#print(g_pf3d7_dict_neighbor.keys())
#print(g_cpar_dict_neighbor)




