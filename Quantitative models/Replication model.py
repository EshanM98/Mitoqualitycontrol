# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 23:47:42 2023

@author: Eshan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functions import *


data = pd.read_excel(r'G:\My Drive\Code\Data\All Experimental Data.xlsx')
deg_model = pd.read_excel(r'G:\My Drive\Code\Data\Degradation model (90 gen 90% f0).xlsx')
dataprev_all = np.concatenate((100*data['Prev'].dropna(), data['Soma']))
datanext_all = np.concatenate((100*data['Next'].dropna(), data['Ovary']))

x = np.linspace(0,100,100)
bins = 20
binthresh = 3
prevgen, nextgen,groups_prev,groups_next = binning(bins,binthresh,dataprev_all,datanext_all)
# prevgen = dataprev_all
# nextgen = datanext_all
mut_gen = np.linspace(0.01,0.9999,1000)
lamb_nocomp= [1.22,1.505,1.646,1.727,1.777,1.808]
lamb_r1 = 1
lamb_r2 = [1,1,1.1616,1.2929,1.3939,1.4646]
mtc_fin = [280,560,1120,2240,4480,8960]

f_comp_r1 = []
f_comp_r2 = []
f_nocomp = []

for i in range(len(mut_gen)):
    f = mut_gen[i]
    f_nextgen_nocomp = rep_nocomp(f,lamb_nocomp[0],mtc_fin[0])
    f_nextgen_comp_r1 = rep_r1(f,lamb_r1,mtc_fin[5])[0]
    f_nextgen_comp_r2 = rep_r2(f,lamb_r2[5],mtc_fin[5])[0]
    
    f_comp_r1.append(f_nextgen_comp_r1)
    f_comp_r2.append(f_nextgen_comp_r2)
    f_nocomp.append(f_nextgen_nocomp)  
    
sse_nocomp = sse_rep_nocomp(lamb_nocomp[0],prevgen,nextgen,mtc_fin[0])
sse_comp_r1 = sse_rep_r1_alt(lamb_r1,prevgen,nextgen,n[5])
sse_comp_r2 = sse_rep_r2(lamb_r2[5],prevgen,nextgen,mtc_fin[5])



fig, ax = plt.subplots(figsize = (15,15))
ax.plot(x,x,'--',label = 'Offspring = parent',ls = (0,(5,10)),lw=4)
ax.scatter(np.asarray(dataprev_all),np.asarray(datanext_all),color ='grey',s=90)
#ax.plot(100*mut_gen,100*np.asarray(f_alt_r1),lw = 6,)
ax.plot(100*mut_gen,100*np.asarray(f_nocomp),label = 'No complementation',color = 'g',lw = 6)
ax.plot(100*mut_gen,100*np.asarray(f_comp_r1),label = 'Region 1 with complementation',color = 'g',lw = 6,ls = '-.')
ax.plot(100*mut_gen,100*np.asarray(f_comp_r2),label = 'Region 2A with complementation',color = 'g',lw = 6, ls = '--')
ax.legend(fontsize = 25,handlelength=3.5)
ax.set_xlabel('Mutant percent in previous generation',fontsize = 34)
ax.set_ylabel('Mutant percent in next generation',fontsize = 34)
ax.tick_params(axis='both',which='major',labelsize=30)
ax.set_xlim(0,100)
ax.set_ylim(0,100)
plt.show()


print('SSE R1 = ' + str(sse_comp_r1))
print('SSE R2a = ' +str(sse_comp_r2))
print('SSE no comp = ' + str(sse_nocomp))
