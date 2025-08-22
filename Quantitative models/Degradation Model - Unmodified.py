# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 23:13:50 2023

@author: Eshan
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functions import *
import matplotlib

data = pd.read_excel(r'G:\My Drive\Code\Data\All Experimental Data.xlsx')
deg_model = pd.read_excel(r'G:\My Drive\Code\Data\Degradation model (90 gen 90% f0).xlsx')
dataprev_all = np.concatenate((100*data['Prev'].dropna(), data['Soma']))
datanext_all = np.concatenate((100*data['Next'].dropna(), data['Ovary']))
ma_dataprev = 100*data['Prev_Ma'].dropna()
ma_datanext = 100*data['Next_Ma'].dropna()
hill_dataprev = 100*data['Prev_Hill'].dropna()
hill_datanext = 100*data['Next_Hill'].dropna()
hurd_dataprev = 100*data['Prev_Hurd'].dropna()
hurd_datanext = 100*data['Next_Hurd'].dropna()
hurd_soma = data['Soma'].dropna() #prev
hurd_ovary = data['Ovary'].dropna() #next


bins = 20
binthresh = 3
prevgen, nextgen, groups_prev,groups_next = binning(bins,binthresh,dataprev_all,datanext_all)


f_reg2A= []
f_reg1 = []
f_int = []
gen = np.linspace(0,0.999,100)
p = 0.323

for i in range(len(gen)):
    f_next = deg2A(gen[i])[0]
    f_reg2A.append(f_next)


for i in range(len(gen)):
    f_next = deg1(gen[i])[0]
    f_reg1.append(f_next)

for i in range(len(gen)):
    f_next = inter_reg(gen[i], p)
    f_int.append(f_next)
x = np.linspace(0,100,100)

sse_r1,sse_r2,sse_int = sse(prevgen,nextgen,p)

fig,ax  = plt.subplots(figsize = (15,15))
# ax.scatter(np.asarray(hill_dataprev),np.asarray(hill_datanext),color = 'indigo',s= 90)
# ax.scatter(np.asarray(ma_dataprev),np.asarray(ma_datanext),color = 'mediumorchid',s= 90)
# ax.scatter(np.asarray(hurd_dataprev),np.asarray(hurd_datanext),color = 'royalblue',s= 90)
# ax.scatter(np.asarray(hurd_soma),np.asarray(hurd_ovary),color = 'pink',s= 90)
ax.plot(x,x,ls = (0,(5,10)), label = 'Offspring = parent',lw = 4)
ax.scatter(np.asarray(dataprev_all),np.asarray(datanext_all),color = 'grey', s = 90)
ax.plot(100*gen,100*np.asarray(f_reg1),color = 'brown',label = 'Region 1',lw = 4)
ax.plot(100*gen,100*np.asarray(f_reg2A),color ='m',label = 'Region 2A',lw = 4)
ax.plot(100*gen,100*np.asarray(f_int),label = 'Interpolated region',lw=6,color = 'r')
ax.set_xlabel('Mutant percent in previous generation',fontsize = 35)
ax.set_ylabel('Mutant percent in next generation',fontsize = 35)
ax.legend(fontsize = 30, handlelength=3.5)
ax.tick_params(axis = 'both',which='major',labelsize = 30)
ax.set_xlim(0,100)
ax.set_ylim(0,100)
plt.show()


print('SSE r1 = ' + str(round(sse_r1,3)))
print('SSE r2a = ' +str(round(sse_r2,3)))
print('SSE int = ' + str(round(sse_int,3)))