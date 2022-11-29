#!/usr/bin/env python
# coding: utf-8

def lim_targets(t):
	lis=[]
	t_unique=t.drop_duplicates('A-genes')
	for x in t_unique['A-genes']:
		tl=t.loc[t['A-genes'] == x]
		if len(tl) > 3:
			lis.append(x)
	filt=pd.DataFrame({'A-genes': lis}).drop_duplicates()
	out=pd.merge(t, filt, on='A-genes', how='inner')

	return out