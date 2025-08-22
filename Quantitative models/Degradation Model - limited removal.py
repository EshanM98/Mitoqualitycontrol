# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 23:46:54 2023

@author: eshan_user
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

mut_gen = np.linspace(0.99,0.01,100)
f_int = [mut_gen[0]]
p= 0.323
for i in range(len(mut_gen)):
    f_next = inter_reg(mut_gen[i], p)
    f_int.append(f_next)

f_intrem = [mut_gen[0]]
p_rem = 0.202
for i in range(len(mut_gen)):
    f_next = deg2A_int(mut_gen[i],p_rem)
    f_intrem.append(f_next)

f_r2a_onerem = [mut_gen[0]]
f_r2a_tworem = [mut_gen[0]]
f_r2a_threerem = [mut_gen[0]]
f_r1_onerem = [mut_gen[0]]
f_r1_tworem = [mut_gen[0]]
f_r1_threerem = [mut_gen[0]]
for i in range(len(mut_gen)):
    f_next = deg2A_onerem(mut_gen[i])
    f_r2a_onerem.append(f_next)
    
    f_next = deg2A_tworem(mut_gen[i])
    f_r2a_tworem.append(f_next)
    
    
    f_next = deg1_onerem(mut_gen[i])
    f_r1_onerem.append(f_next)
    
    f_next = deg1_tworem(mut_gen[i])
    f_r1_tworem.append(f_next)
    
    f_next = deg1_threerem(mut_gen[i])
    f_r1_threerem.append(f_next)
    
bins = 20
binthresh = 3
prevgen, nextgen, groups_prev,groups_next = binning(bins,binthresh,dataprev_all,datanext_all)
# prevgen = dataprev_all
# nextgen = datanext_all

sse_r1,sse_r2,sse_int = sse(prevgen,nextgen,p)

x = np.linspace(0,100,100)

prev_r1 = np.asarray(deg_model['Region 1 prediction previous gen'])
next_r1 = np.asarray(deg_model['Region 1 prediction next gen'])


prev_r2 = np.asarray(deg_model['Region 2A prediction previous gen'])
next_r2 = np.asarray(deg_model['Region 2A prediction next gen'])


prev_int = np.asarray(f_int[:-1])
next_int = np.asarray(f_int[1:])


prev_r2_one = np.asarray(f_r2a_onerem[:-1])
next_r2_one = np.asarray(f_r2a_onerem[1:])


prev_r2_two = np.asarray(f_r2a_tworem[:-1])
next_r2_two = np.asarray(f_r2a_tworem[1:])



prev_r1_one = np.asarray(f_r1_onerem[:-1])
next_r1_one = np.asarray(f_r1_onerem[1:])


prev_r1_two = np.asarray(f_r1_tworem[:-1])
next_r1_two = np.asarray(f_r1_tworem[1:])


prev_r1_three = np.asarray(f_r1_threerem[:-1])
next_r1_three = np.asarray(f_r1_threerem[1:])



sse_2a_onerem = []
sse_2a_tworem = []
sse_1_onerem = []
sse_1_tworem = []
sse_2a_int = []

for i in range(len(prevgen)):
    y0 = deg2A_onerem(prevgen[i]/100)
    y = nextgen[i]/100
    err = ((y-y0)/(prevgen[i]/100))**2
    sse_2a_onerem.append(err)
    
    y0 = deg2A_tworem(prevgen[i]/100)
    y = nextgen[i]/100
    err = ((y-y0)/(prevgen[i]/100))**2
    sse_2a_tworem.append(err)
    
    y0 = deg1_onerem(prevgen[i]/100)
    y = nextgen[i]/100
    err = ((y-y0)/(prevgen[i]/100))**2
    sse_1_onerem.append(err)
    
    y0 = deg1_tworem(prevgen[i]/100)
    y = nextgen[i]/100
    err = ((y-y0)/(prevgen[i]/100))**2
    sse_1_tworem.append(err)
    
    y0 = deg2A_int(prevgen[i]/100,p_rem)
    y = nextgen[i]/100 
    err = ((y-y0)/(prevgen[i]/100))**2
    sse_2a_int.append(err)
    


fig,ax = plt.subplots(figsize = (15,15))
cmap = matplotlib.cm.get_cmap('inferno')
ax.plot(x,x,lw = 4,ls = (0,(5,10)),label = 'Offspring = parent')
ax.scatter(np.asarray(dataprev_all),np.asarray(datanext_all),color = 'grey',s=90)
#ax.scatter(np.asarray(hill_dataprev),np.asarray(hill_datanext),color = 'indigo',s= 90)
#ax.scatter(np.asarray(ma_dataprev),np.asarray(ma_datanext),color = 'mediumorchid',s= 90)
#ax.scatter(np.asarray(hurd_dataprev),np.asarray(hurd_datanext),color = 'royalblue',s= 90)
#ax.scatter(np.asarray(hurd_soma),np.asarray(hurd_ovary),color = 'pink',s= 90)
#ax.plot(100*prev_r1, 100*next_r1,'--',  label ='Region 1. SSE = ' + str(sse_r1), color = 'brown',lw = 4)
#ax.plot(100*prev_r2,100*next_r2,'--',  label = 'Region 2A. SSE = ' + str(sse_r2),color = 'm',lw = 4)
# ax.plot(100*mut_gen,100*next_int,'--',label = 'Interpolated degradation model. SSE = ' + str(sse_int))
ax.plot(100*mut_gen,100*next_r1_one,lw = 4,color = 'brown',ls = '--',label = ' Region 1 one nucleoid removal.')
ax.plot(100*mut_gen,100*next_r1_two,lw=4,color = 'brown',label = ' Region 1 one/two nucleoid removal.')
#ax.plot(100*mut_gen,100*next_r1_three,linewidth = 4,color = 'brown',ls='-.',label = 'Region 1 three/two/one nucleoid removed.')
ax.plot(100*mut_gen,100*next_r2_one,linewidth = 4,color = 'm',ls = '--',label = 'Region 2a one nucleoid removed.')
ax.plot(100*mut_gen,100*next_r2_two,linewidth = 4,color = 'm',label = 'Region 2a one/two nucleoid removed. ')
#ax.plot(100*mut_gen,100*np.asarray(f_intrem[1:]), linewidth = 6,color = 'r',label ='Interpolated Region 2a degradation.')
ax.set_xlabel('Mutant percent in previous generation', fontsize = 35)
ax.set_ylabel('Mutant percent in next generation', fontsize = 35)
ax.legend(fontsize = 23,handlelength=3.4)
ax.tick_params(axis='both',which='major',labelsize = 30)
ax.set_xlim(0,100)
ax.set_ylim(0,100)
plt.show()

# print('SSE r1 = ' + str(sse_r1))
# print('SSE r2a = ' +str(sse_r2))
# print('SSE int = ' + str(sse_int))
print('SSE r1_onerem = ' + str(round(sum(sse_1_onerem),2)))
print('SSE r1_tworem = ' + str(round(sum(sse_1_tworem),2)))
print('SSE r2a_onerem = ' +str(round(sum(sse_2a_onerem),2)))
print('SSE r2a_tworem = ' +str(round(sum(sse_2a_tworem),2)))
print('SSE int = ' + str(round(sum(sse_2a_int),2)))