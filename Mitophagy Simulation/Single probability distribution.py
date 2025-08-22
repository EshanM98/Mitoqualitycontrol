# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 15:31:42 2024

@author: eshan_user
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom
import pandas as pd
from matplotlib.lines import Line2D
from func_gillespie_stochastic import *
from prob_func import *
from scipy import integrate
import time
start_time = time.time()
## Initialize size distribution
# a1 = height of peak
#mu =  center of peak
sigma = 0.25 # st_dev
size = 100000
nuc_prob = [0.739,0.226,0.035]
nuc_count = [1,2,3]
x= np.linspace(0.75,5,size)
f_lim = 0.75

probs1,size1,y_norm1 = init(x,nuc_prob[0],nuc_count[0],sigma,f_lim)
probs2,size2,y_norm2 = init(x,nuc_prob[1],nuc_count[1],sigma,f_lim)
probs3,size3,y_norm3 = init(x,nuc_prob[2],nuc_count[2],sigma,f_lim)




size_dist = [size1,size2,size3]
prob_dist = [probs1,probs2,probs3]

##Initialize simulation
protein_init = 10000 #Max protein count
protein_split = 1 #if 0 - wt only prod, if 1- mutant only prod, if 0.5 - both prod
prodrate = protein_init/19  #production of protein per day
tau = 2 #half life
back_mito = 0.01 # cannot be greater than 1/19
mito = 1 #if 1, active mitophagy occurs
mtc_count = 100
end_time = 20
threshold = 0.25
mut_init = 0.5 #Starting quantity of mutant
wt_count = protein_init*(1-mut_init) #Starting quantity of wildtype proteins
end = 0.95
active_prod = 0
drop = 0.5
bite_size = 0.25
std_bite = 0.02

#Arrays to track data
probs_track = []
t_track = []
t_exp_track  = []
kdp_avg_rec= []
true_bites_rec = []
bite_draws = []

fig,ax = plt.subplots(figsize = (15,15))
colors = ['r','b','g']
for j in range(len(size_dist)):
    mtcm_wcount,mtcm_mcount,mtcm_totcount,tm_tot,t_fin,mitorate,decay,bite_size_rec,bite_draw_rec = gillespie_stoch(protein_init,mut_init,prodrate, mtc_count,end_time,tau,protein_split,mito,back_mito,bite_size,threshold,active_prod,drop,size_dist[j],prob_dist[j],std_bite)
    bite_draws.append(bite_draw_rec)
    true_bites_rec.append(bite_size_rec)
    kdp_arr = np.asarray(decay) + np.asarray(mitorate)*np.asarray(bite_draw_rec)
    kdp = sum(kdp_arr)/len(kdp_arr)
    kdp_avg_rec.append(kdp)
    t_exp = -np.log(threshold/(1-mut_init))/(kdp) #0 prodrate expectation time
    t_exp_track.append(t_exp)
        
        
    probs = []
       
    for i in range(len(t_fin)):
        probs.append((i+1)/len(t_fin))
    probs_track.append(probs)
        
    t_fin.sort()
    t_track.append(t_fin)
        
        
    #idx1_fin = probs_track[0].index(end) #When nuc 1 reaches 95%
    #t1_end = t_track[0][idx1_fin]
    colors = ['r','b','g']
    ax.plot(t_fin,probs,label = str(nuc_count[j]),color = colors[j], lw = 3)
    #ax.plot([t_exp,t_exp],[0,1],ls = (0,(5,10)),color = colors[j],lw = 2)
    

t_all = []
for i in range(len(t_track)):
    t_all.extend(t_track[i])
t_all.sort()



p_all = []
p = 0
for j in range(len(t_all)):
    for i in range(len(t_track)):
        if t_all[j] in t_track[i]:
            idx = t_track[i].index(t_all[j])
            #p = probs_track[i][idx]*nuc_probs[i]
            if probs_track[i][idx] == 0:
                p = p
            else:
                p = p + 1/(len(t_fin))*nuc_prob[i]
            p_all.append(p)
            
bins = 100
bin_width = round(len(p_all)/bins)

p_binned = [sum(p_all[i:i+bin_width])/bin_width for i in range(0,len(p_all),bin_width)]
t_binned = [sum(t_all[i:i+bin_width])/bin_width for i in range(0,len(t_all),bin_width)]
    
kdp_avg = 0
for i in range(len(kdp_avg_rec)):
    kdp_avg = kdp_avg + kdp_avg_rec[i]*nuc_prob[i]
t_exp = -np.log(threshold/(1-mut_init))/(kdp_avg)


#ax.plot(t_binned,p_binned,color = 'm', lw = 5)
#ax.plot([t_exp,t_exp],[0,1],ls = (0,(5,10)),color = 'm', lw = 5)
idx1_fin = probs_track[0].index(end) #When nuc 1 reaches 95%
t1_end = t_track[0][idx1_fin]
ax.plot([t1_end,t1_end],[0,1],ls = (0,(5,10)), color = 'salmon', lw= 5)
ax.set_xlabel('Time (days)',fontsize = 30)
ax.set_ylabel('Threshold Probability',fontsize = 30)
ax.legend(title = 'Nucleoid count',title_fontsize = 30, fontsize = 26)
ax.plot()
#ax.set_title('St dev: ' +str(sigma) + ' Bite size: ' + str(bite_size) + ' Bite size st dev: ' + str(std_bite),fontsize = 30)
ax.tick_params(axis='both',which='major',labelsize=30)
ax.set_ylim(bottom = 0)
ax.set_xlim([0,10])
plt.show()
