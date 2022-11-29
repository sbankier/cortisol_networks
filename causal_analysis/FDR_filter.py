#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np

#import findr networks
res=pd.read_csv('_pijs_gassist_GR_full_A-it_cis-anchor_alt_1Mb', delimiter='\t')

#filter networks to selected FDR threshold
def filt_net(net, trd, test):
    steps=np.arange(0.5, 1, 0.01)
    fdr_l=[]
    for s in steps:
        fdr=1-net.loc[net[test] >= s][test].mean()
        fdr_l.append(fdr)
    fdr_filt=pd.DataFrame({test: steps, 'FDR': fdr_l})
    sel=fdr_filt.iloc[(fdr_filt['FDR']-trd).abs().argsort()[:1]].values[0][0]
    selr=np.round(sel, 2)
    net_fdr=net.loc[net[test] >= selr].sort_values(test, ascending=False)
    
    return net_fdr

fdr10=filt_net(res_filt, 0.1, 'p2p5')
