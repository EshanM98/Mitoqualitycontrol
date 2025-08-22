# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 13:59:24 2023

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
prevgen, nextgen,groups_prev,groups_next = binning(bins,binthresh,dataprev_all,datanext_all)


lamb = np.linspace(1,2,100)
mtc_fin = []
num_of_reps = np.linspace(1,6,20)
for i in range(len(num_of_reps)):
    fin = 140*(2**num_of_reps[i])
    mtc_fin.append(fin)

#No Complementation Error minimization
lamb_min_nocomp = []
sse_lamb_nocomp = []
for k in range(len(mtc_fin)):
    sse_sum = []
    for i in range(len(lamb)):
        sse_nocomp = sse_rep_nocomp(lamb[i],prevgen,nextgen,mtc_fin[k])
        sse_sum.append(sse_nocomp)
    lamb_min = lamb[sse_sum.index(min(sse_sum))]
    lamb_min_nocomp.append(lamb_min)
    sse_lamb_nocomp.append(min(sse_sum))


#Region 1 Error Minimization
lamb_min_r1 = []
sse_lamb_r1 = []
for k in range(len(mtc_fin)):
    sse_sum = []
    for i in range(len(lamb)):
        sse_comp = sse_rep_r1(lamb[i], prevgen, nextgen, mtc_fin[k])[1]
        sse_sum.append(sse_comp)
    lamb_min = lamb[sse_sum.index(min(sse_sum))]
    lamb_min_r1.append(lamb_min)
    sse_lamb_r1.append(min(sse_sum))


#Region 2 Error minimization
lamb_min_r2 = []
sse_lamb_r2 = []
for k in range(len(mtc_fin)):
    sse_sum = []
    for i in range(len(lamb)):
        sse_comp = sse_rep_r2(lamb[i], prevgen, nextgen, mtc_fin[k])[1]
        sse_sum.append(sse_comp)
    lamb_min = lamb[sse_sum.index(min(sse_sum))]
    lamb_min_r2.append(lamb_min)
    sse_lamb_r2.append(min(sse_sum))


fig,ax = plt.subplots(figsize = (15,15))
ax.plot(num_of_reps,lamb_min_r2,label = 'Region 2a',color = 'g',ls = '--', lw = 5)
ax.plot(num_of_reps,lamb_min_r1,label = 'Region 1',color = 'g',ls = '-.',lw = 5)
ax.plot(num_of_reps,lamb_min_nocomp,label = 'No Complementation',lw = 3, color = 'g')
ax.set_xlabel('Replication Cycles',fontsize = 34)
ax.set_ylabel('Lambda Value',fontsize =34)
ax.tick_params(axis='both',which='major',labelsize=30)
ax.legend(fontsize=30,handlelength=3.5)
ax.set_ylim(0.9,2)
ax.set_xlim(1,6)
plt.show()

fig, ax = plt.subplots(figsize = (15,15))
ax.plot(num_of_reps,sse_lamb_r2,label = 'Region 2a',color = 'g',ls = '--',lw = 5)
ax.plot(num_of_reps,sse_lamb_r1,label = 'Region 1',color = 'g', ls = '-.', lw = 5)
ax.plot(num_of_reps,sse_lamb_nocomp,label = 'No Complementation',lw = 3,color = 'g')
ax.set_xlabel('Replication Cycles',fontsize = 34)
ax.set_ylabel('Sum of Squared Errors',fontsize = 34)
ax.set_xlim(1,6)
ax.tick_params(axis='both',which='major',labelsize=30)
ax.legend(fontsize=30,handlelength=3.5)
plt.show()






