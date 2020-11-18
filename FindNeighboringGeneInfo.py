from pybedtools import BedTool
import os
import re
from cor_distance import g_cor_dict
import random

#read in corresponding gff files
g_output_cor_relation="/Users/yiranli/Dropbox/lncRNA_article/data/cor/cor_relationship.tab"
os.chdir("/Users/yiranli/Dropbox/lncRNA_article/data/cor")
os.system('''cat /Users/yiranli/Dropbox/lncRNA_article/lncRNA_checked_backed.gff|awk '{if($3=="exon") print}'|sort -k1,1 -k4,4n|uniq>tmp_final_noncoding_transcript.gff''')
os.system('''cat /Users/yiranli/Desktop/CryptoDB-43_CparvumIowaII.gff|awk '{if($3=="mRNA") print}'|sort -k1,1 -k4,4n|uniq>tmp_mRNA_transcript.gff''')

#build gene relationship database by bedtools
lncrna_db=BedTool("tmp_final_noncoding_transcript.gff")
mrna_db=BedTool("tmp_mRNA_transcript.gff")
upstream=lncrna_db.closest(mrna_db,io=True,id=True,D="a").saveas("tmp_upstream.tab")
downstream=lncrna_db.closest(mrna_db,io=True,iu=True,D="a").saveas("tmp_downstream.tab")
antisense=lncrna_db.intersect(mrna_db,wao=True,S=True).saveas("tmp_antisense.tab")
os.system('''cat tmp_upstream.tab|awk -F "\t" '{if($3=="exon" && $19<-1) print $9"\t"$18"\t"$19}'|sort|uniq >tmp_upstream.index''') 
os.system('''cat tmp_downstream.tab|awk -F "\t" '{if($3="exon" && $19!=-1) print $9"\t"$18"\t"$19}'|sort|uniq >tmp_downstream.index''') 
os.system('''cat tmp_antisense.tab|awk -F "\t" '{if($3=="exon" && $19!=0) print $9"\t"$18"\t"$19}'|sort|uniq >tmp_antisense.index''')

g_lncRNA_id=[]
for line in lncrna_db:
	line=str(line)
	lncRNA_id=re.search(r'Parent=.*',line).group().replace("Parent=","").replace("\tID=","")+"-RA"
	g_lncRNA_id.append(lncRNA_id)

# for each lncRNA, find the sense mRNA, upstream and downstream mRNA information
dict_upstream={}
dict_downstream={}
dict_antisense={}

def make_relation_dict(input_file):
	with open(input_file) as f:
		a_dict={}
		for line in f:
			try:
				key=re.search(r'Parent=.*\tID=',line).group().replace("Parent=","").replace("\tID=","")+"-RA"
				val=re.search(r'[Cc]gd\d_\d+-RA',line).group().replace("Cgd","cgd")
				a_dict[key]=val
			except:
				pass			
	return a_dict

dict_upstream=make_relation_dict("tmp_upstream.index")
dict_downstream=make_relation_dict("tmp_downstream.index")
dict_antisense=make_relation_dict("tmp_antisense.index")

#output lncRNA and the neighboring mRNA information table
with open(g_output_cor_relation,"w") as write_to:
	write_to.write("relationship\tcor\n")
	for gene in g_lncRNA_id:
		try:
			write_to.write('lncRNA_upstream_mRNA\t{}\n'.format(g_cor_dict[(gene,"mRNA_"+dict_upstream[gene])]))
		except:
			pass
		try:
			write_to.write('lncRNA_downstream_mRNA\t{}\n'.format(g_cor_dict[(gene,"mRNA_"+dict_downstream[gene])]))
		except:
			pass
		try:
			write_to.write('lncRNA_sense_mRNA\t{}\n'.format(g_cor_dict[(gene,"mRNA_"+dict_antisense[gene])]))
		except:
			pass
	gene_random_gene_cor=random.sample(list(g_cor_dict.values()),k=10000)
	for i in gene_random_gene_cor:
		write_to.write('gene_random_gene\t{}\n'.format(i))










	
