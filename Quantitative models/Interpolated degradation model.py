# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 15:50:17 2023

@author: eshan_user
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functions import *

data = pd.read_excel(r'G:\My Drive\Code\Data\All Experimental Data.xlsx')
deg_model = pd.read_excel(r'G:\My Drive\Code\Data\Degradation model (90 gen 90% f0).xlsx')
dataprev_all = np.concatenate((100*data['Prev'].dropna(), data['Soma']))
datanext_all = np.concatenate((100*data['Next'].dropna(), data['Ovary']))

bins = 20
binthresh = 3
prevgen, nextgen, groups_prev,groups_next = binning(bins,binthresh,dataprev_all,datanext_all)
# prevgen = dataprev_all
# nextgen = datanext_all

p = np.linspace(0,1,100)
sse_2a_int = []
sse_int = []
for i in range(len(p)):
    sse_sumint = sse(prevgen,nextgen,p[i])[2]
    sse_int.append(sse_sumint)

for i in range(len(p)):
    sse_2a = []
    for j in range(len(prevgen)):
        y0 = deg2A_int(prevgen[j]/100,p[i])
        y = nextgen[j]/100 
        err = ((y-y0)/(prevgen[j]/100))**2
        sse_2a.append(err)
    sse = sum(sse_2a)
    sse_2a_int.append(sse)
        
fig,ax = plt.subplots(figsize = (15,15))
ax.plot(p,sse_2a_int, label = 'Interpolation of selective degradation model in region 2a')
ax.plot(p,sse_int,label = 'Interpolation of unmodified degradation model')
ax.set_xlabel('P-value')
ax.set_ylabel('SSE')
ax.legend()
ax.set_title('Changes in SSE between regions of degradation model')
plt.show()

print('Size-selective degradation p-val: ' + str(p[sse_2a.index(min(sse_2a))]))
print('Regular degradation p-val: ' + str(p[sse_int.index(min(sse_int))]))
