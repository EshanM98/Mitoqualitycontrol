# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 11:47:22 2024

@author: eshan_user
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom
import pandas as pd
from matplotlib.lines import Line2D
from func_gillespie_stochastic import *
from prob_func import *
import seaborn as sns
import time
start_time = time.time()


#Set up size distribution

sigma = 0.25
size = 100000
nuc_prob = [0.739,0.211,0.05]
nuc_count = [1,2,3]
x = np.linspace(0.1,5,size)
f_lim = 0.75


probs1,size1,y_norm1 = init(x,nuc_prob[0],nuc_count[0],sigma,f_lim)
probs2,size2,y_norm2 = init(x,nuc_prob[1],nuc_count[1],sigma,f_lim)
probs3,size3,y_norm3 = init(x,nuc_prob[2],nuc_count[2],sigma,f_lim)

#rates - lambda of gillespie algorithm

tau = 2 #half life

back_mito = 0.01
mito = 1


bites = np.linspace(0.05,0.5,10)
std_bite = 0.02
active_prod = 0
#active_prod = np.linspace(0,1,10)
#drop = 1/2
drop = np.linspace(0.1,0.75,10)
size_dist = [size1,size2,size3]
prob_dist = [probs1,probs2,probs3]

mtc_count = 1000
protein_init = 10000 #Max protein count, start at 50%
prodrate = protein_init/19
end_time = 20
threshold = 0.25
mut_init = 0.5


end = 0.95

protein_split = 1 

p2_stop_prob = []
p3_stop_prob = []
prod_step_track = 0
for m in range(len(drop)):
    p2_stop_rec = []
    p3_stop_rec = []

    bite_track = 0
    for k in range(len(bites)):
        
        
        mtcm_wcount,mtcm_mcount,mtcm_totcount,tm_tot,t_fin1,mitorate,decay,bite_size_true,bite_draw_rec = gillespie_stoch(protein_init,mut_init,prodrate, mtc_count,end_time,tau,protein_split,mito,back_mito,bites[k],threshold,active_prod,drop[m],size_dist[0],prob_dist[0],std_bite)
            
        t_fin1.sort()
            
        probs1 = []
            
        for i in range(len(t_fin1)):
            probs1.append(i/len(t_fin1))
            
        idx1_fin = probs1.index(end) #When nuc 1 reaches 95%
        t1_end = t_fin1[idx1_fin]
            
        mtcm_wcount,mtcm_mcount,mtcm_totcount,tm_tot,t_fin2,mitorate,decay,bite_size_true,bite_draw_rec = gillespie_stoch(protein_init,mut_init,prodrate, mtc_count,t1_end,tau,protein_split,mito,back_mito,bites[k],threshold,active_prod,drop[m],size_dist[1],prob_dist[1],std_bite)
        probs2 = []
        for i in range(len(t_fin2)):
            probs2.append(i/mtc_count)
        if len(probs2)>0:
            p2_stop_rec.append(probs2[-1])
        else:
            p2_stop_rec.append(0)
        
        mtcm_wcount,mtcm_mcount,mtcm_totcount,tm_tot,t_fin3,mitorate,decay,bite_size_true,bite_draw_rec = gillespie_stoch(protein_init,mut_init,prodrate, mtc_count,t1_end,tau,protein_split,mito,back_mito,bites[k],threshold,active_prod,drop[m],size_dist[2],prob_dist[2],std_bite)
        probs3 = []
        for i in range(len(t_fin3)):
            probs3.append(i/mtc_count)
        if len(probs3)>0:
            p3_stop_rec.append(probs3[-1])
        else:
            p3_stop_rec.append(0)
            
        bite_track +=1  
        print('Current position: Active Production = ' +str(active_prod[m]) +'Bite size = ' + str(bites[k]))
    p2_stop_prob.append(p2_stop_rec)
    p3_stop_prob.append(p3_stop_rec)
    
    prod_step_track += 1
    print("Active Production step is:", prod_step_track)
    
p2_df = pd.DataFrame(data = p2_stop_prob[::-1], columns = np.round(bites,2), index = np.round(active_prod[::-1],2))
s2 = sns.heatmap(p2_df,cbar_kws={'label': 'Threshold Probability'})
#s2.set(xlabel = 'Bite Size',ylabel = 'Fraction of Production Rate')
s2.set(xlabel = 'Bite Size',ylabel = 'Fraction of proteins remaining')
#plt.title('Probability of mito 2 when 95% of mito 1 reach threshold')
plt.title('Nucleoid count: 2')
plt.show()


p3_df = pd.DataFrame(data = p3_stop_prob[::-1],columns = np.round(bites,2), index = np.round(active_prod[::-1],2))
s3=sns.heatmap(p3_df,cbar_kws={'label': 'Threshold Probability'})
#s3.set(xlabel = 'Bite Size',ylabel = 'Fraction of Production Rate')
s3.set(xlabel = 'Bite Size',ylabel = 'Fraction of proteins remaining')
#plt.title('Probability of mito 3 when 95% of mito 1 reach threshold')
plt.title('Nucleoid count: 3')
plt.show()
    
    
print("Code runtime:", time.time() - start_time)

