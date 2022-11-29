import pandas as pd
import numpy as np
import math
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

#import adjusted METSIM data
sf_met_adj=pd.read_csv('GSE135134_METSIM_subcutaneousAdipose_RNAseq_log-TPM1_n434_PCA1_adj-PC1.txt', delimiter='\t')
sf_met_adj.rename(columns={'Gene.stable.ID':'Gene stable ID'}, inplace=True)

def pca_plt(ct, pca1, pca2, save, title):

    #apply standard scaler to dataset features
    ct_x = ct.drop(['Gene stable ID'], 1).dropna().values
    ct_xt = StandardScaler().fit_transform(ct_x)

    #reduce data across 2 principal components
    ct_pca = PCA(n_components=2)
    ct_principalComponents = ct_pca.fit_transform(ct_x)
    ct_principalDf = pd.DataFrame(data = ct_principalComponents
                 , columns = ['principal component 1', 'principal component 2'])

    #remove outliers from PCA
    if pca1 == 0:
        pass
    else:
        ct_principalDf=ct_principalDf.sort_values('principal component 1', ascending=False).iloc[pca1:]

    if pca1 == 0:
        pass
    else:
        ct_principalDf=ct_principalDf.sort_values('principal component 2', ascending=False).iloc[pca2:]

    #plot PCA and save output
    sns.lmplot(x='principal component 1', y='principal component 2', fit_reg=False,
                 palette='mako_r', data=ct_principalDf.dropna(), height=7)
    plt.tight_layout()
    if save == True:
        plt.savefig(title+'.png')
    plt.show()

pca_plt(sf_metf, 0, 0, True, 'GSE135134_METSIM_subcutaneousAdipose_RNAseq_log-tmp1_n434')
pca_plt(sf_met_adj, 0, 0, True, 'GSE135134_METSIM_subcutaneousAdipose_RNAseq_log-tmp1_adj-PC1_n434')

def get_pca(ct):

    #apply standard scaler to dataset features
    ct_x = ct.drop(['Gene stable ID'], 1).dropna().values
    ct_xt = StandardScaler().fit_transform(ct_x)

    #reduce data across 2 principal components
    ct_pca = PCA(n_components=2)
    ct_principalComponents = ct_pca.fit_transform(ct_x)
    ct_principalDf = pd.DataFrame(data = ct_principalComponents
                 , columns = ['principal component 1', 'principal component 2'])

    return ct_principalDf

sf_pca=get_pca(sf_metf)
sf_metf['principal component 1']=sf_pca['principal component 1'].to_list()
cols = sf_metf.columns.tolist()
cols = cols[-1:] + cols[:-1]
sf_pca_out=sf_metf[cols]