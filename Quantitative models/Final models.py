
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
dataprev_all = np.concatenate((100*data['Prev'].dropna(), data['Soma']))
datanext_all = np.concatenate((100*data['Next'].dropna(), data['Ovary']))


x = np.linspace(0,100,100)
bins = 20
binthresh = 3
prevgen, nextgen,groups_prev,groups_next = binning(bins,binthresh,dataprev_all,datanext_all)

#replication initiate
mut_gen = np.linspace(0.01,0.9999,100)
lamb_nocomp= [1.22,1.556,1.768,1.879]
lamb_r1 = 1
lamb_r2 = [1,1,1.2,1.35,1.4,1.4646]
mtc_fin = [280,560,1120,2240,4480,8960]

f_comp_r1 = []
f_comp_r2 = []
f_nocomp = []
f_alt_r1 = []
for i in range(len(mut_gen)):
    f = mut_gen[i]
    f_nextgen_nocomp = rep_r2(f,lamb_nocomp[0],mtc_fin[0])[0]
    f_nextgen_comp_r1 = rep_r1(f,lamb_r1,mtc_fin[5])[1]
    f_nextgen_comp_r2 = rep_r2(f,lamb_r2[5],mtc_fin[5])[1]
    
    
    f_comp_r1.append(f_nextgen_comp_r1)
    f_comp_r2.append(f_nextgen_comp_r2)
    f_nocomp.append(f_nextgen_nocomp)  
    
    
sse_nocomp = sse_rep_r2(lamb_nocomp[0],prevgen,nextgen,mtc_fin[0])[0]
sse_comp_r1 = sse_rep_r1(lamb_r1,prevgen,nextgen,mtc_fin[5])[1]
sse_comp_r2 = sse_rep_r2(lamb_r2[5],prevgen,nextgen,mtc_fin[5])[1]


print('SSE R1 = ' + str(sse_comp_r1))
print('SSE R2a = ' +str(sse_comp_r2))
print('SSE no comp = ' + str(sse_nocomp))

#Degradation initiate
mut_gen = np.linspace(0.01,0.9999,100)

f_intrem = [mut_gen[0]]
p_rem = 0.202
for i in range(len(mut_gen)):
    f_next = deg2A_int(mut_gen[i],p_rem)
    f_intrem.append(f_next)



x = np.linspace(0,100,100)



p = 0.75

sse_2a_onerem = []
sse_2a_tworem = []
sse_1_onerem = []
sse_1_tworem = []
sse_2a_int = []

for i in range(len(prevgen)):
    y0 = deg2A_int(prevgen[i]/100,p_rem)
    y = nextgen[i]/100 
    err = ((y-y0)/(prevgen[i]/100))**2
    sse_2a_int.append(err)



fig, ax = plt.subplots(figsize = (15,15))
ax.plot(x,x,'--',label = 'Offspring = parent',ls = (0,(5,10)),lw=4)
ax.scatter(np.asarray(dataprev_all),np.asarray(datanext_all),color = 'grey',s = 90)
ax.plot(100*mut_gen,100*np.asarray(f_intrem[1:]),color = 'r',label ='Interpolated degradation model',linewidth =6)
ax.plot(100*mut_gen,100*np.asarray(f_comp_r2),color = 'g',label = 'Region 2A replication model',linewidth = 6,ls = '--')
ax.set_xlabel('Mutant percent in previous generation',fontsize = 34)
ax.set_ylabel('Mutant percent in next generation', fontsize = 34)
ax.set_xlim(0,100)
ax.set_ylim(0,100)
ax.legend(fontsize = 25,handlelength=3.5)
ax.tick_params(axis='both',which='major',labelsize = 30)
plt.show()