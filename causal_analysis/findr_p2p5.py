#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import findr
import numpy as np

""""Iterates aross list of genes with best cis-eQTL and calculates all 5 tests from Findr 
using best eQTLs as causal anchors. Manual calculation of pij_gassist, pij_gassist_trad 
and P2*P5 test combination and returns 3 DataFrames as pandas panel and exports output as pickle"""

#state specific tissue e.g. LIV, SKLM, MAM, AOR, VAF, SF, BLOOD and file containing genes and eQTLs of interest
t='VAF'
genes=pd.read_csv(t+'_GR_target_priority_cis-eQTLs_1Mb', delimiter='\t')
genes=genes.loc[genes['Findr_score'] >= 0.75]

#import expression data
t_mrna=pd.read_csv(t+'.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN', delimiter=' ')
p=t_mrna['id']
id=[]
for x in p:
	id.append(x[0:15])
t_mrna=t_mrna.drop('id', axis=1)
t_mrna['Gene stable ID']=id
cols=t_mrna.columns.tolist()
cols=cols[-1:] + cols[:-1]
t_mrna=t_mrna[cols]

trad_out=[]
novel_out=[]
alt_out=[]
for index, row in genes.iterrows():
	#import gene annotation and get eQTLs
	nm=pd.read_csv(t+'/'+t+'_gene_annotation', delimiter='\t')
	eQTLs=pd.DataFrame({'rs_number': [row[0]], 'Gene_name': [row[1]], 'Findr_score': [row[2]], 'Gene stable ID': [row[3]]})

	#filter to get eQTL genotypes
	df_list=[]
	eQTLs.rename(columns={'rs_number':'marker_id'}, inplace=True)
	rs=eQTLs[['marker_id']]
	for cnk in pd.read_csv(t+'/'+t+'.STARNET.v3.genotype.dose.MAF-0.05',
			sep='\t', iterator=True, chunksize=200000):
			mer=pd.merge(cnk, rs, on='marker_id', how='inner')
			df_list.append(mer)
	trans_geno=pd.concat(df_list)

	#obtain expression data from gene list
	trans_nm=pd.DataFrame(eQTLs['Gene stable ID'])
	trans_mrna=pd.merge(trans_nm, t_mrna, on='Gene stable ID', how='left')

	#remove duplicate genes so there is expression data for specificed eQTLs
	mrna_index=t_mrna.loc[t_mrna['Gene stable ID'].isin(trans_mrna['Gene stable ID'].tolist())]
	mrnai=mrna_index.reset_index()
	mrna_index_sorted=pd.merge(trans_nm, mrnai, on='Gene stable ID', how='left')
	mrna_index_v=mrna_index_sorted['index'].values
	t_mrna_sorted=pd.concat([t_mrna.iloc[mrna_index_v,:],
		t_mrna.drop(mrna_index_v, axis=0)], axis=0)

	#remove columns and convert to format for findr
	dg=trans_geno.drop('marker_id', axis=1)
	d=trans_mrna.drop('Gene stable ID', axis=1)
	dt=t_mrna_sorted.drop('Gene stable ID', axis=1)

	dgi=dg.values
	di=d.values
	dti=dt.values

	dg_input=dgi.astype(np.uint8)
	d_input=di.astype(np.float32)
	dt_input=dti.astype(np.float32)

	#run novel causal inference test in Findr
	l=findr.lib()
	ans=l.pijs_gassist(dg_input,d_input,dt_input,na=None,nodiag=True)
	ans

	#create x and y labels
	xlab=t_mrna_sorted['Gene stable ID']
	ylab=pd.merge(trans_mrna, genes, on='Gene stable ID', how='left').drop_duplicates(subset='Gene_name')

	#generate and orientate dataframes
	p1=pd.DataFrame(ans['p1'])
	p1t=p1.T.reset_index()
	p2=pd.DataFrame(ans['p2'], columns=xlab)
	p2t=p2.T.reset_index()
	p3=pd.DataFrame(ans['p3'], columns=xlab)
	p3t=p3.T.reset_index()
	p4=pd.DataFrame(ans['p4'], columns=xlab)
	p4t=p4.T.reset_index()
	p5=pd.DataFrame(ans['p5'], columns=xlab)
	p5t=p5.T.reset_index()

	#replace Ensembl ID with genes names
	p1t.columns=['Gene stable ID']+list(ylab['Gene_name'])
	p2t.columns=['Gene stable ID']+list(ylab['Gene_name'])
	p3t.columns=['Gene stable ID']+list(ylab['Gene_name'])
	p4t.columns=['Gene stable ID']+list(ylab['Gene_name'])
	p5t.columns=['Gene stable ID']+list(ylab['Gene_name'])

	p2_comp=pd.merge(p2t, nm, on='Gene stable ID', how='inner')
	p2_comp=p2_comp.drop(['Unnamed: 0', 'Transcript stable ID', 'Chromosome/scaffold name',
	 'Gene start (bp)', 'Gene end (bp)', 'Gene stable ID'], axis=1)
	p3_comp=pd.merge(p3t, nm, on='Gene stable ID', how='inner')
	p3_comp=p3_comp.drop(['Unnamed: 0', 'Transcript stable ID', 'Chromosome/scaffold name',
	 'Gene start (bp)', 'Gene end (bp)', 'Gene stable ID'], axis=1)
	p4_comp=pd.merge(p4t, nm, on='Gene stable ID', how='inner')
	p4_comp=p4_comp.drop(['Unnamed: 0', 'Transcript stable ID', 'Chromosome/scaffold name',
	 'Gene start (bp)', 'Gene end (bp)', 'Gene stable ID'], axis=1)
	p5_comp=pd.merge(p5t, nm, on='Gene stable ID', how='inner')
	p5_comp=p5_comp.drop(['Unnamed: 0', 'Transcript stable ID', 'Chromosome/scaffold name',
	'Gene start (bp)', 'Gene end (bp)', 'Gene stable ID'], axis=1)

	#reorder columns so Gene Name comes first
	cols=p2_comp.columns.tolist()
	cols=cols[-1:] + cols[:-1]
	p2_comp=p2_comp[cols]
	p3_comp=p3_comp[cols]
	p4_comp=p4_comp[cols]
	p5_comp=p5_comp[cols]

	#remove gene columns for test calculation
	p2_m=p2_comp.drop('Gene_name', axis=1)
	p3_m=p3_comp.drop('Gene_name', axis=1)
	p4_m=p4_comp.drop('Gene_name', axis=1)
	p5_m=p5_comp.drop('Gene_name', axis=1)

	#Get results of P2*P5 test
	p2p5=p2_m*p5_m
	p2p5['Gene_name']=p2_comp['Gene_name']
	p2p5=p2p5[cols]
	p2p5=p2p5.set_index('Gene_name')
	p2p5=p2p5.stack().reset_index()
	p2p5.columns=['B-genes', 'A-genes', 'Findr_score']
	p2p5_df=p2p5[['A-genes', 'B-genes', 'Findr_score']]
	alt_out.append(p2p5_df)

	#Get results for pij_gassist_trad
	trad=p2_m*p3_m
	trad['Gene_name']=p2_comp['Gene_name']
	trad=trad[cols]
	trad=trad.set_index('Gene_name')
	trad=trad.stack().reset_index()
	trad.columns=['B-genes', 'A-genes', 'Findr_score']
	trad_df=trad[['A-genes', 'B-genes', 'Findr_score']]
	trad_out.append(trad_df)

	#Get results for pij_gassist
	gas=(0.5*((p2_m*p5_m)+p4_m))
	gas['Gene_name']=p2_comp['Gene_name']
	gas=gas[cols]
	gas=gas.set_index('Gene_name')
	gas=gas.stack().reset_index()
	gas.columns=['B-genes', 'A-genes', 'Findr_score']
	gas_df=gas[['A-genes', 'B-genes', 'Findr_score']]
	novel_out.append(gas_df)

#export data as dataframes
alt_exp=pd.concat(alt_out)
trad_exp=pd.concat(trad_out)
novel_exp=pd.concat(novel_out)
p2p5_exp.to_csv(t+'_pijs_gassist_GR_full_A-it_cis-anchor_alt_1Mb', sep='\t', index=None)
trad_exp.to_csv(t+'_pijs_gassist_GR_full_A-it_cis-anchor_trad_1Mb', sep='\t', index=None)
novel_exp.to_csv(t+'_pijs_gassist_GR_full_A-it_cis-anchor_novel_1Mb', sep='\t', index=None)
