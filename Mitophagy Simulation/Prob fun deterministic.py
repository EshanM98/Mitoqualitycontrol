# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 09:58:41 2023

@author: eshan_user
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.special as scspec
from prob_func import *




f_size = [1,2,3]
size_prob = [0.739,0.226,0.035]
probs = [0.739,0.211,0.05]
bite = 0.25
max_t = 70
P_fin = np.zeros(max_t+1).tolist()
k_rec = 0
k_track = []
fig, ax = plt.subplots(figsize = (15,15))

for i in range(len(f_size)):
    f_min = f_size[0]
    f_stop = 0.25
    f_remove = bite* f_min/f_size[i]
    print(f_remove)
    fissionrate = 1
    t = np.linspace(0,max_t,max_t+1)
    if f_remove < 1:
        k_start = np.ceil(np.log(f_stop)/(np.log(1-f_remove)))
    
    print(k_start)
    P_rec = []
    k_track.append(k_start)
    for j in range(len(t)):
        lamb = fissionrate * t[j]
        if lamb < 40:
            p = 0
            k_max = 10*k_start
            for k in range(int(k_start),int(k_max)):
                if k < 160:
                                
                
                    p_new = poisson(k,lamb)
                else:
                    p_new = norm(k,lamb)              
                            
                p = p + p_new
                            
            P_rec.append(p)
                        
        else:
            k_start = np.floor(np.log(f_stop)/np.log(1-f_remove))
            p_new = 1 - (0.5 + 0.5*scspec.erf((k_start-lamb)/np.sqrt(2*lamb)))
            P_rec.append(p_new)
                        
            
    color = ['r','b','g','m','c','k']
            
    P_fin = np.asarray(P_fin) + np.asarray(P_rec) * size_prob[i]
    k_rec = k_rec + k_start*size_prob[i]
    
    ax.plot(t,P_rec, label = str(f_size[i]), color = color[i],lw = 3)
    ax.plot([k_start/fissionrate,k_start/fissionrate],[0,1],'--', color = color[i],lw = 2)



ax.plot([(k_rec/fissionrate),(k_rec/fissionrate)],[0,1],linestyle ='dashed',color = 'm',lw = 5)
ax.plot(t,P_fin, linewidth = 4, color = 'm')
ax.set_xlabel('Time',fontsize = 30)
ax.set_ylabel('Threshold probability',fontsize = 30)
ax.legend(title = 'Nucleoid count',title_fontsize = 30, fontsize = 26)
ax.plot()
#ax.set_title('St dev: ' +str(sigma) + ' Bite size: ' + str(bite_size) + ' Bite size st dev: ' + str(std_bite),fontsize = 30)
ax.tick_params(axis='both',which='major',labelsize=30)
ax.set_ylim(bottom = 0)
ax.set_xlim([0,25])
plt.show()





