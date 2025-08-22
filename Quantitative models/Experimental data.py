# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:24:56 2024

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


ma_fly1 = ma_dataprev.loc[0:9].to_list()
ma_fly1.append(ma_datanext.iloc[9])
ma_fly2 = ma_dataprev.loc[10:19].to_list()
ma_fly2.append(ma_datanext.iloc[19])
ma_fly3 = ma_dataprev.loc[20:29].to_list()
ma_fly3.append(ma_datanext.iloc[29])
ma_gen = np.arange(0,(len(ma_fly1)),1)

ma_fly1 = pd.DataFrame({'gen': ma_gen,'ma_fly' : ma_fly1})
ma_fly2 = pd.DataFrame({'gen':ma_gen, 'ma_fly' : ma_fly2})
ma_fly3 = pd.DataFrame({'gen':ma_gen, 'ma_fly' : ma_fly3})
ma_gen_all = pd.concat([ma_fly1,ma_fly2,ma_fly3])


hill_fly = hill_dataprev.loc[0::].to_list()
hill_fly.append(hill_datanext.iloc[-1])
hill_gen = np.arange(0,len(hill_fly),1)

hill_dat = {'gen': hill_gen, 'hill_fly' : hill_fly}
hill_gen_all = pd.DataFrame(data = hill_dat)

hurd_fly1 = hurd_dataprev.loc[0:3].to_list()
hurd_fly1.append(hurd_datanext.iloc[3])
hurd_fly2 = hurd_dataprev.loc[4:6].to_list()
hurd_fly2.append(hurd_datanext.iloc[6])
hurd_fly3 = hurd_dataprev.loc[7:8].to_list()
hurd_fly3.append(hurd_datanext.iloc[8])
hurd_fly4 = [hurd_dataprev.iloc[9]]
hurd_fly4.append(hurd_datanext.iloc[9])
hurd_fly5 = hurd_dataprev.loc[10:12].to_list()
hurd_fly5.append(hurd_datanext.iloc[12])
hurd_fly6 = hurd_dataprev.loc[13:15].to_list()
hurd_fly6.append(hurd_datanext.iloc[15])
hurd_fly7 = hurd_dataprev.loc[16:18].to_list()
hurd_fly7.append(hurd_datanext.iloc[18])
hurd_fly8 = hurd_dataprev.loc[19:21].to_list()
hurd_fly8.append(hurd_datanext.iloc[21])
hurd_fly9 = hurd_dataprev.loc[22:23].to_list()
hurd_fly9.append(hurd_datanext.iloc[23])
hurd_fly10 = hurd_dataprev.loc[24:27].to_list()
hurd_fly10.append(hurd_datanext.iloc[27])
hurd_fly11 = hurd_dataprev.loc[28:29].to_list()
hurd_fly11.append(hurd_datanext.iloc[29])

hurd_dat1 = pd.DataFrame({'gen':np.arange(0,len(hurd_fly1),1),'hurd_fly':pd.Series(hurd_fly1)}) 
hurd_dat2 =pd.DataFrame({'gen':np.arange(0,len(hurd_fly2),1),'hurd_fly':pd.Series(hurd_fly2)})
hurd_dat3 = pd.DataFrame({'gen':np.arange(0,len(hurd_fly3),1),'hurd_fly':pd.Series(hurd_fly3)})
hurd_dat4 = pd.DataFrame({'gen':np.arange(0,len(hurd_fly4),1),'hurd_fly':pd.Series(hurd_fly4)})
hurd_dat5 = pd.DataFrame({'gen':np.arange(0,len(hurd_fly5),1), 'hurd_fly':pd.Series(hurd_fly5)})
hurd_dat6 = pd.DataFrame({'gen':np.arange(0,len(hurd_fly6),1),'hurd_fly':pd.Series(hurd_fly6)})
hurd_dat7 = pd.DataFrame({'gen':np.arange(0,len(hurd_fly7),1),'hurd_fly':pd.Series(hurd_fly7)})
hurd_dat8 = pd.DataFrame({'gen':np.arange(0,len(hurd_fly8),1),'hurd_fly':pd.Series(hurd_fly8)})
hurd_dat9 = pd.DataFrame({'gen':np.arange(0,len(hurd_fly9),1),'hurd_fly':pd.Series(hurd_fly9)})
hurd_dat10 = pd.DataFrame({'gen':np.arange(0,len(hurd_fly10),1),'hurd_fly':pd.Series(hurd_fly10)})
hurd_dat11 = pd.DataFrame({'gen':np.arange(0,len(hurd_fly11),1),'hurd_fly':pd.Series(hurd_fly11)})
hurd_gen_all = pd.concat([hurd_dat1,hurd_dat2,hurd_dat3,hurd_dat4,hurd_dat5,hurd_dat6,hurd_dat7,hurd_dat8,hurd_dat9,hurd_dat10,hurd_dat11])


gen_soma = [0]*len(hurd_soma)
gen_ovary = [1]*len(hurd_ovary)

hurd_gen_soma = pd.DataFrame({'gen':gen_soma,'fly':hurd_soma})
hurd_gen_ovary = pd.DataFrame({'gen':gen_ovary,'fly':hurd_ovary})

hurd_so_gen = []
hurd_so_fly = []
for i in range(len(hurd_gen_soma)):
    hurd_so_gen.append(hurd_gen_soma.iloc[i,0])
    hurd_so_fly.append(hurd_gen_soma.iloc[i,1])
    
    hurd_so_gen.append(hurd_gen_ovary.iloc[i,0])
    hurd_so_fly.append(hurd_gen_ovary.iloc[i,1])

hurd_gen_so_all = pd.DataFrame({'gen':hurd_so_gen,'fly':hurd_so_fly})

figure,ax = plt.subplots(figsize = (15,15))
ax.scatter(hurd_gen_so_all['gen'],hurd_gen_so_all['fly'],color = 'pink',s=90)
ax.scatter(ma_gen_all['gen'],ma_gen_all['ma_fly'],color = 'mediumorchid',s=90)
ax.scatter(hill_gen_all['gen'],hill_gen_all['hill_fly'],color = 'indigo',s=90)
ax.scatter(hurd_gen_all['gen'],hurd_gen_all['hurd_fly'],color = 'royalblue',s=90)
ax.set_xlabel('Generation',fontsize = 35)
ax.set_ylabel('Percentage of mutant genome',fontsize = 35)
ax.set_ylim(0,100)
#ax.legend(fontsize = 23,handlelength=3.4)
ax.tick_params(axis='both',which='major',labelsize = 30)
plt.show()


figure, ax = plt.subplots(figsize = (15,15))
ax.scatter(np.asarray(hill_dataprev),np.asarray(hill_datanext),color = 'indigo',s= 90)
ax.scatter(np.asarray(ma_dataprev),np.asarray(ma_datanext),color = 'mediumorchid',s= 90)
ax.scatter(np.asarray(hurd_dataprev),np.asarray(hurd_datanext),color = 'royalblue',s= 90)
ax.scatter(np.asarray(hurd_soma),np.asarray(hurd_ovary),color = 'pink',s= 90)
ax.set_xlabel('Mutant percent in previous generation',fontsize = 35)
ax.set_ylabel('Mutant percent in next generation',fontsize = 35)
ax.set_ylim(0,100)
#ax.legend(fontsize = 23,handlelength=3.4)
ax.tick_params(axis='both',which='major',labelsize = 30)
plt.show()