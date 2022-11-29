#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from scipy import stats
import math
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import distance
from scipy.cluster import hierarchy

#general functions
def flip_df(df):
    df=df.T.reset_index()
    df.columns = df.iloc[0]
    df=df.iloc[1:]
    return df

def swap_cols(df):
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    return df


# # Correlation heatmaps
# 
# This analysis examines correlations between network targets of GR-regulated gene networks using gene expression data from STAGE and METSIM. Heatmaps show correlations between network targets using STARNET vs METSIM or STAGE. Genes have been clustered according to STARNET gene expression data and the same clustering has been applied to the corresponding replication datasets.

#import gene networks
sf=pd.read_csv('SF_pijs_gassist_alt_cis-anchor_0.75_1Mb_CORNET-GWAMA_trans-genes_15FDR_GR-targets.tsv', delimiter='\t')
liv=pd.read_csv('LIV_pijs_gassist_alt_cis-anchor_0.75_1Mb_CORNET-GWAMA_trans-genes_15FDR_GR-targets.tsv', delimiter='\t')
vaf=pd.read_csv('VAF_pijs_gassist_alt_cis-anchor_0.75_1Mb_CORNET-GWAMA_trans-genes_15FDR_GR-targets.tsv', delimiter='\t')

sf10=sf.loc[sf['Findr_score'] >= 0.87]
liv10=liv.loc[liv['Findr_score'] >= 0.855]
vaf10=vaf.loc[vaf['Findr_score'] >= 0.86]

#get absolute correlation values as matrix
def get_corr(gene, net, exp):
    net_f=net.loc[net['A-genes'] == gene]
    net_met=pd.merge(net_f[['B-genes']], exp, left_on='B-genes', right_on='Gene_name').drop_duplicates('Gene_name')
    sf_t=flip_df(net_met.drop(['Gene stable ID', 'B-genes'], axis=1))
    corr=sf_t.astype(float).drop('Gene_name', axis=1).corr().abs().reset_index().fillna(0)
    corr.rename(columns={0:'Gene_name'}, inplace=True)

    return corr


ir10_met=get_corr('IRF2', sf10, sf_met_adj)

#import STARNET network correlations
ir10_star_ori=pd.read_csv('SF_STARNET_IRF2_10FDR_net_correlations.tsv', sep='\t')

def process_star_corr(ori):
    ori.rename(columns={'0':'Gene_name'}, inplace=True)
    return ori

ir10_star=process_star_corr(ir10_star_ori)

#import STAGE SF data
nm_ref=pd.read_csv('STAGE/NM_ref', delimiter='\t')

def process_stage_exp(tis, ref):
    stage_full=pd.read_csv('STAGE_'+tis+'_RMA.txt', delimiter='\t')
    stage_full.rename(columns={'Unnamed: 0':'NM_number'}, inplace=True)
    stage_anno=pd.merge(stage_full, nm_ref, on='NM_number', how='inner').drop('mrna', axis=1)
    stage_mer=pd.merge(stage_anno, ref, on='Gene_name', how='inner').drop_duplicates('Gene_name')
    stage=swap_cols(stage_mer)
    stage.rename(columns={'NM_number':'Gene stable ID'}, inplace=True)
    return stage

sf_stage=process_stage_exp('SubcutFat', sf_ref)
ir10_stage=get_corr('IRF2', sf10, sf_stage)

def filt_corr(df_big, df_small):
    filt_rows=pd.merge(df_small[['Gene_name']], df_big, on='Gene_name', how='inner')
    cols=df_small.columns
    filt_cols=filt_rows[cols]
    
    return filt_cols

ir10_star_stage=filt_corr(ir10_star, ir10_stage)

def apply_clustering(star_df, df):
    
    correlations_array = np.asarray(star_df.set_index('Gene_name'))

    row_linkage = hierarchy.linkage(
        distance.pdist(correlations_array), method='average')

    col_linkage = hierarchy.linkage(
        distance.pdist(correlations_array.T), method='average')
    
    star_ref=pd.DataFrame(star_df['Gene_name']).reset_index()
    new_ids=pd.DataFrame(hierarchy.leaves_list(row_linkage), columns=['index'])
    new_index=pd.merge(new_ids, star_ref, on='index', how='inner')
    re_rows=pd.merge(new_index[['Gene_name']], df, on='Gene_name', how='inner')
    re_cols=re_rows[re_rows['Gene_name'].tolist()]
    
    return re_cols

#IRF2 SF
ir10_star_clust=apply_clustering(ir10_star, ir10_star)
ir10_met_clust=apply_clustering(ir10_star, ir10_met)
ir10_stage_clust=apply_clustering(ir10_star_stage, ir10_stage)
ir10_star_stage_clust=apply_clustering(ir10_star_stage, ir10_star_stage)

def plot_heatmap(corr1, corr2, lab1, lab2, save, title):

    # Generate a mask for the upper triangle
    mask1 = np.triu(np.ones_like(corr1, dtype=bool))
    mask2 = np.triu(np.ones_like(corr2, dtype=bool)).T

    # Set up the matplotlib figure
    plt.figure(figsize=(20,15))
    plt.rcParams.update({'font.size': 26})

    # Draw the heatmap with the mask
    ax=sns.heatmap(corr1, square=True, cmap="YlGnBu", mask=mask1, 
                   cbar_kws={'label': '\n'+lab1+' Pearson absolute value'}, vmin=0, vmax=1)
    ax.tick_params(labelbottom=False, labelleft=False)
    ax.tick_params(axis='both', which='both', length=0)
   # plt.title(title)


    sns.heatmap(corr2, square=True, cmap="magma_r", mask=mask2, 
                cbar_kws={'label': '\n'+lab2+' Pearson absolute value'}, vmin=0, vmax=1)
    ax.tick_params(labelbottom=False, labelleft=False) 
    ax.tick_params(axis='both', which='both', length=0)

    if save == True:
        plt.savefig('correlation_heatmaps/'+title+'_heatmap.png')

plot_heatmap(ir10_met_clust, ir10_star_clust, 'METSIM', 'STARNET', True, 'SF_IRF2_corr_STARNET_METSIM')

# # Obtain correlations between random and network targets
# 
# Comparing distribution of correlations between network targets to distribution of random gene set of the same size across different datasets.

sf_star_ori=pd.read_csv('SF.exp.mat.F.gene_filt.EDAseq.gc.lm.or.RQN', delimiter='\t')
star_ref=pd.read_csv('gene_annotation_GRCh37_biomart', delimiter='\t')
sf_star=pd.merge(star_ref[['Gene stable ID', 'Gene name']].drop_duplicates(), sf_star_ori, on='Gene stable ID', how='inner')
sf_star.rename(columns={'Gene name':'Gene_name'}, inplace=True)

#get absolute correlation values as matrix from random sample
def get_corr_rand(gene, net, exp):
    net_f=net.loc[net['A-genes'] == gene]
    rn=len(net_f)
    net_exp=exp.drop_duplicates('Gene_name').sample(n=rn, random_state=100)
    net_t=flip_df(net_exp.drop(['Gene stable ID'], axis=1))
    corr=net_t.astype(float).drop('Gene_name', axis=1).corr().abs().reset_index().fillna(0)
    corr.rename(columns={0:'Gene_name'}, inplace=True)

    return corr

def rmv_diag(mat):
    matnan=mat.replace(1.0,np.NaN).drop('Gene_name', axis=1).values
    mat_nodiag=matnan[np.logical_not(np.isnan(matnan))]
    
    return mat_nodiag

def make_df(d1, d2, d1l, d2l, d3l):
    d1_df=pd.DataFrame({'corr': d1})
    d2_df=pd.DataFrame({'corr': d2})
    d1_df['label']=[d1l]*len(d1_df)
    d2_df['label']=[d2l]*len(d2_df)
    d3_df=pd.concat([d1_df, d2_df])
    d3_df['data']=[d3l]*len(d3_df)
    
    
    return d3_df

def get_kw(t, gene, target):
    tar=t.loc[t['label'] == target]['corr'].values
    ran=t.loc[t['label'] == 'random']['corr'].values
    kw=stats.kruskal(tar, ran)
    t['kw_stat']=[kw[0]]*len(t)
    t['p-value']=[kw[1]]*len(t)
    
    return t

#IRF2 SF
ir10_met_rand=get_corr_rand('IRF2', sf10, sf_met_adj)
ir10_stage_rand=get_corr_rand('IRF2', sf10, sf_stage)
ir10_star_rand=get_corr_rand('IRF2', sf10, sf_star)

ir10_met_df=make_df(rmv_diag(ir10_met), rmv_diag(ir10_met_rand), '10FDR', 'random', 'METSIM')
ir10_stage_df=make_df(rmv_diag(ir10_stage), rmv_diag(ir10_stage_rand), '10FDR', 'random', 'STAGE')
ir10_star_df=make_df(rmv_diag(ir10_star), rmv_diag(ir10_star_rand), '10FDR', 'random', 'STARNET')

ir10_met_kw=get_kw(ir10_met_df, 'IRF2', '10FDR')
ir10_stage_kw=get_kw(ir10_stage_df, 'IRF2', '10FDR')
ir10_star_kw=get_kw(ir10_star_df, 'IRF2', '10FDR')

ir10_box=pd.concat([ir10_met_kw, ir10_stage_kw, ir10_star_kw])

def plt_box(plot, gene, tar, save):
    kw=plot.drop_duplicates('data').sort_values('data', ascending=False)['p-value'].values
    l=[]
    for x in kw:
        if x == 0:
            kwps='< 1.000e-300'
        else:
            kwps=f'{x:.3e}'
        l.append(kwps)
    plt.figure(figsize=(12,10))
    ax = sns.boxplot(data=plot, x="data", y="corr", hue='label')
    plt.ylabel('Pearson absolute value', fontsize=22)
    plt.xlabel('')
    plt.yticks(fontsize=22)
    plt.xticks(fontsize=22)
    plt.legend(fontsize=20)
    lab=plot.drop_duplicates('data').sort_values('data', ascending=False)['data'].values
    if len(l) == 3:
        ax.set_xticklabels([lab[0]+'\n\n(p = '+str(l[0])+')', lab[1]+'\n\n(p = '+str(l[1])+')', lab[2]+'\n\n(p = '+str(l[2])+')'])
    else:
        ax.set_xticklabels([lab[0]+'\n\n(p = '+str(l[0])+')', lab[1]+'\n\n(p = '+str(l[1])+')'])
    if save == True:
        plt.savefig(gene+'_'+tar+'_correlation_boxplot.png', bbox_inches='tight')
    plt.show()


plt_box(ir10_box, 'IRF2', '10FDR', True)
