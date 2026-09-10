import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import pandas as pd, numpy as np, os
W='wgs/04_ancestry'; out=open(f'{W}/summary.txt','w')
for tag,kg,dy in [('GLOBAL','kg.proj.sscore','target.proj.sscore'),('EAS','eas.proj.sscore','eas.target.proj.sscore')]:
    k=pd.read_csv(f'{W}/{kg}',sep='\t').rename(columns={'#IID':'IID'}); d=pd.read_csv(f'{W}/{dy}',sep='\t')
    pcs=[c for c in k.columns if c.endswith('_AVG')][:4]; D=d[pcs].values[0]
    cen=k.groupby('Population')[pcs].mean(); dist=np.sqrt(((cen-D)**2).sum(1)).sort_values()
    k['dist']=np.sqrt(((k[pcs].values-D)**2).sum(1))
    s=f'== {tag} (PC1-4): {NAME_EN} {np.round(D,4)}\n{dist.head(6).round(4).to_string()}\nnearest 15 individuals: {k.nsmallest(15,"dist").Population.value_counts().to_dict()}\n'
    # {NAME_EN}'s percentile position within CHB/CHS along PC1/PC2 of the EAS PCA
    if tag=='EAS':
        for p in ['CHB','CHS','JPT','KHV','CDX']:
            sub=k[k.Population==p]; s+=f'{p}: PC1 pct {(sub[pcs[0]]<D[0]).mean()*100:.0f}, PC2 pct {(sub[pcs[1]]<D[1]).mean()*100:.0f}\n'
    print(s); out.write(s)
