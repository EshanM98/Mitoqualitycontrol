# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 19:45:56 2024

@author: eshan_user
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom
import pandas as pd
from matplotlib.lines import Line2D
from func_gillespie_stochastic import *
from prob_func import *
import time
start_time = time.time()


#Set up size distribution

sigma = 0.5
size = 100000
nuc_prob = [0.739,0.211,0.05]
nuc_count = [1,2,3]
x = np.linspace(0.1,5,size)
f_lim = 0.75

prob1,size1,y_int1 = init(x,nuc_prob[0],nuc_count[0],sigma,f_lim)
prob2,size2,y_int2 = init(x,nuc_prob[1],nuc_count[1],sigma,f_lim)
prob3,size3,y_int3 = init(x,nuc_prob[2],nuc_count[2],sigma,f_lim)


#rates - lambda of gillespie algorithm

tau = 2 #half life

back_mito = 0.01
mito = 1
bites = np.linspace(0.05,0.5,10)

size_dist = [size1,size2,size3]
prob_dist = [prob1,prob2,prob3]

mtc_count = 1000
protein_init = 10000 #Max protein count, start at 50%
prodrate = protein_init/19
end_time = 20
threshold = 0.25
mut_init = 0.5
stdev_bite = 0.02

end = 0.95
drop = 1/2
protein_split = 1 
active_prod = 0

p2_stop_rec = []
p3_stop_rec = []
bite_track = 0
fig,ax = plt.subplots(figsize = (15,15))
for k in range(len(bites)):
    mtcm_wcount,mtcm_mcount,mtcm_totcount,tm_tot,t_fin1,mitorate,decay,bite_size_true,bite_draw_rec1 = gillespie_stoch(protein_init,mut_init,prodrate, mtc_count,end_time,tau,protein_split,mito,back_mito,bites[k],threshold,active_prod,drop,size_dist[0],prob_dist[0],stdev_bite)
    
    t_fin1.sort()
    
    probs1 = []
    
    for i in range(len(t_fin1)):
        probs1.append(i/len(t_fin1))
    
    idx1_fin = probs1.index(end) #When nuc 1 reaches 95%
    t1_end = t_fin1[idx1_fin]
    
    mtcm_wcount,mtcm_mcount,mtcm_totcount,tm_tot,t_fin2,mitorate,decay,bite_size_true,bite_draw_rec2 = gillespie_stoch(protein_init,mut_init,prodrate, mtc_count,t1_end,tau,protein_split,mito,back_mito,bites[k],threshold,active_prod,drop,size_dist[1],prob_dist[1],stdev_bite)
    probs2 = []
    
    for i in range(len(t_fin2)):
        probs2.append(i/mtc_count)
    if len(probs2)>0:
        p2_stop_rec.append(probs2[-1])
    else:
        p2_stop_rec.append(0)
    
    mtcm_wcount,mtcm_mcount,mtcm_totcount,tm_tot,t_fin3,mitorate,decay,bite_size_true,bite_draw_rec3 = gillespie_stoch(protein_init,mut_init,prodrate, mtc_count,t1_end,tau,protein_split,mito,back_mito,bites[k],threshold,active_prod,drop,size_dist[2],prob_dist[2],stdev_bite)
    probs3 = []
    
    for i in range(len(t_fin3)):
        probs3.append(i/mtc_count)
    if len(probs3)>0:
        p3_stop_rec.append(probs3[-1])
    else:
        p3_stop_rec.append(0)
        
    bite_track += 1
    print("Bite step is: " + str(bite_track) + '. Bite Size is: ' + str(bites[k]))
    
    
    
    
ax.plot([0,1],[0.95,0.95],label = '1',color = 'r',lw = 5)
ax.plot(bites,p2_stop_rec,label = '2',color = 'b', lw = 5)
ax.plot(bites,p3_stop_rec,label = '3',color = 'g',lw=5)
ax.plot([0.25,0.25],[0,1],ls = (0,(5,10)), color = 'salmon', lw= 4)
ax.set_xlabel('Bite size',fontsize = 30)
ax.set_ylabel('Threshold Probability',fontsize = 30)
ax.legend(title = 'Nucleoid count',title_fontsize = 30,fontsize = 26)
#ax.set_title('st dev:' +str(sigma),fontsize = 30)
#ax.set_title('st dev of bite: ' + str(stdev_bite),fontsize = 30)
ax.tick_params(axis='both',which='major',labelsize=30)
ax.set_xlim([0,0.5])
ax.set_ylim([0,1])
plt.show()

