#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import findr
import numpy as np

"""Script takes as input STARNET tissue type runs secondary linkage test
in Findr to examine association between given SNPs and genes within the tissue. Returns
Findr scores as a matrix of rs numbers against gene names and value indicating association
i.e 1=strong association 0=independent"""

#specify tissue e.g. LIV, SKLM, MAM, AOR, VAF, SF, Blood
t='LIV'

#import expression data
t_mrna=pd.read_csv(+t+'.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN', delimiter=' ')
p=t_mrna['id']
id=[]
for x in p:
	id.append(x[0:15])
t_mrna=t_mrna.drop('id', axis=1)
t_mrna['Gene stable ID']=id
cols=t_mrna.columns.tolist()
cols=cols[-1:] + cols[:-1]
t_mrna=t_mrna[cols]

#import genotype data and gene annotation
t_geno=pd.read_csv(t+'.STARNET_genotype.MAF-0.05_SERPINA6_100Kb', delimiter='\t')
t_geno.rename(columns={'0':'rs_number'}, inplace=True)
nm=pd.read_csv(t+'/'+t+'_gene_annotation', delimiter='\t')

#prep d with randomly permutated expression data
geno_shape=np.shape(t_geno)
rand_mrna=np.random.randn(geno_shape[0], geno_shape[1]-1)
d=pd.DataFrame(rand_mrna)

#remove columns and convert to format for findr
dg=t_geno.drop('rs_number', axis=1)
dt=t_mrna.drop('Gene stable ID', axis=1)

dgi=dg.values
di=d.values
dti=dt.values

dg_input=dgi.astype(np.uint8)
d_input=di.astype(np.float32)
dt_input=dti.astype(np.float32)

#calculate all 5 individual tests in Findr
l=findr.lib()
ans=l.pijs_gassist(dg_input,d_input,dt_input,na=None,nodiag=False)
ans

#create x and y labels
xlab=t_mrna['Gene stable ID']
y=t_geno['rs_number']
ylab=pd.DataFrame({'rs_number': y})

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

#Add rs_numbers
p1t.columns=['Gene stable ID']+list(ylab['rs_number'])
p2t.columns=['Gene stable ID']+list(ylab['rs_number'])
p3t.columns=['Gene stable ID']+list(ylab['rs_number'])
p4t.columns=['Gene stable ID']+list(ylab['rs_number'])
p5t.columns=['Gene stable ID']+list(ylab['rs_number'])

#replace ensemble IDs with gene names
ens_ref=nm[['Gene stable ID', 'Gene_name']]
p2_comp=pd.merge(p2t, ens_ref, on='Gene stable ID', how='left')
p3_comp=pd.merge(p3t, ens_ref, on='Gene stable ID', how='left')
p4_comp=pd.merge(p4t, ens_ref, on='Gene stable ID', how='left')
p5_comp=pd.merge(p5t, ens_ref, on='Gene stable ID', how='left')

#reorder columns so Gene Name comes first
cols=p2_comp.columns.tolist()
cols=cols[-1:] + cols[:-1]
p2_comp=p2_comp[cols]
p3_comp=p3_comp[cols]
p4_comp=p4_comp[cols]
p5_comp=p5_comp[cols]

#convert from matrix to table
p2_comp=p2_comp.drop('Gene stable ID', axis=1)
p2_comp=p2_comp.set_index('Gene_name')
p2_comp=p2_comp.stack().reset_index()
p2_comp.columns=['Gene_name', 'rs_number', 'Findr_score']
p2_df=p2_comp[['rs_number', 'Gene_name', 'Findr_score']]

#export output
out.to_csv(t+'_SERPINA6_100Kb_FindrP2_df', sep='\t', index=None)