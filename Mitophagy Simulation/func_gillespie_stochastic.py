# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 00:28:51 2024

@author: eshan_user
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom
import pandas as pd
from matplotlib.lines import Line2D
from prob_func import *




def mito_func(kp,kd_p,f,t,N0):
    Nt = N0*np.exp(-kd_p*t) + (kp/(kd_p)) * (1 - np.exp(-kd_p*t))
    km_active = (kd_p - 1/19)/f
    return Nt, km_active

def steady(km_back,f):
    kd = 1/19 - km_back*f
    return kd


#protein_init - initial protein count
#prodrate - protein production rate (per day)
#mtc_count - number of mitochondria
#end_time - Run time for simulation
#tau - Time (in days) for protein drop to occur
#protein_split - mtDNA presence: 0 means only wtDNA prod, 1 means only mutant prod
#mito - Active mitophagy variable: 0 means active mitophagy is off, 1 means its on
#back_mito - background mitophagy rate
#bite_size - Size of piecemeal mitophagy
#threshold - Minimum wildtype protein fraction at which simulation ends (Fraction of current wildtype/Initial protein count)
#active_prod - fraction of prodrate during active mitophagy
#drop - factor by which protein_init drops in tau
#size_dist - array of mitochondrial sizes
#prob_dist - array of probabilities corresponding to mitochondrial sizes
#bite_stdev - standard deviation of bite size distribution

def gillespie_stoch(protein_init,mut_init, prodrate, mtc_count,end_time,tau,protein_split,mito,back_mito,bite_size,threshold,active_prod,drop,size_dist,prob_dist,bite_stdev):        
    #DataFrames to track protein count over all mutant splits
    mtcm_wcount = pd.DataFrame()
    mtcm_mcount = pd.DataFrame()
    mtcm_totcount = pd.DataFrame()
    tm_tot = pd.DataFrame()
    
    bite_size_rec = []
    bite_draw_rec = []
    decay_rec = []
    mito_rec = []
    kdp_rec = []
    t0 = 0 #start time
    t = t0
    std_bite = bite_stdev #std_dev of gaussian dist of bite sizes
    t_fin = []
    N0 = protein_init   #Protein count
    
    mito_count = 0
    for i in range(mtc_count):
        bite = np.random.normal(bite_size,std_bite)
        #while bite < (bite_size - 0.05) or bite > (bite_size + 0.05) or bite <= 0 or bite >= 0.75:
        while bite <= 0 or bite >= 0.75:
            bite = np.random.normal(bite_size,std_bite)
        
        
        
        
        
        decay = steady(back_mito,bite) #New steady state function (with 1/19 decay timeline)
        decayrate = decay
        decay_rec.append(decayrate)
        #print('decayrate: ' + str(decayrate))

        if mito == 1:
            kp = active_prod*prodrate
            kdp_try = 0.1
            kdp_int = 0.01
            while kdp_int > 0.0001:
                Nt,km_active = mito_func(kp,kdp_try,bite,tau,N0)
                if Nt > N0*drop:
                    kdp_try = kdp_try + kdp_int
                else:
                    kdp_try = kdp_try - kdp_int
                    kdp_int = kdp_int/10
                    kdp_try = kdp_try + kdp_int;
              
            mitorate = km_active + back_mito
            
            #print('kdp: ' +str(kdp_try))
            #print('km_active: '+ str(km_active))
            #print('mitorate: ' + str(mitorate))
         
        
        else:
            kp = prodrate
            mitorate = back_mito
            print(mitorate)
        
        mito_rec.append(mitorate)
         
        
        mtc_size = choose_size(size_dist,prob_dist)
        bite_draw_true = bite/mtc_size
        bite_draw_rec.append(bite_draw_true)
        
        #Starting protein count arrays
        mtc_m_protcount = np.array([int(mut_init * protein_init)]) #initial mutant protein count per mitochondria
        mtc_w_protcount = np.array([int(protein_init - mtc_m_protcount[0])]) #initial wildtype protein count
        
        mtc_tot_protcount = np.array(np.add(mtc_w_protcount,mtc_m_protcount)) #array to track protein count per timestep 
        
        #Arrays to track entire protein count over splits
        wcount = mtc_w_protcount
        mcount = mtc_m_protcount
        totcount = mtc_tot_protcount
        
        t = 0
        t_rec = [0] #array to track time
        bite_sizes = []
    
        while t_rec[-1] < end_time :
            
           bite_new = np.random.normal(bite,std_bite)
           
           while bite_new <= 0 or bite_new >= 0.75:
               bite_new = np.random.normal(bite,std_bite)
           
           bite_size_true = bite_new/mtc_size #fractional piecemeal mitophagy
           bite_sizes.append(bite_size_true)
           
           totalrate = kp  + decayrate*totcount[-1] + mitorate# + splitrate
           delta_split = (1/totalrate)*np.log(1/np.random.rand()) # determine waiting time 
           t = t + delta_split #increment waiting time
           
           prob_prod = kp/totalrate #production probability
           #print(prob_prod)
           prob_decay = decayrate*totcount[-1]/totalrate
           #print(prob_decay)
           prob_mito = mitorate/totalrate
           #print(prob_mito)
           #prob_split = splitrate/totalrate
           #print(prob_split)             
           rand_prob = np.random.rand()
           if rand_prob <= prob_prod:
               rand_prob_2 = np.random.rand()
               if rand_prob_2 <= protein_split:
                    mtc_m_protcount = np.append(mtc_m_protcount,np.add(mtc_m_protcount[-1], 1)) #increase mutant protein count of wildtype mtc by 1
                    mtc_tot_protcount = np.append(mtc_tot_protcount,np.add(mtc_m_protcount[-1], mtc_w_protcount[-1]))#append new total protein value
                    t_rec.append(t)
                    
               else: 
                   mtc_w_protcount = np.append(mtc_w_protcount,np.add(mtc_w_protcount[-1],1))
                   mtc_tot_protcount = np.append(mtc_tot_protcount,np.add(mtc_w_protcount[-1],mtc_m_protcount[-1]))
                   t_rec.append(t)
                   
           elif prob_prod < rand_prob <= prob_prod + prob_decay:
               rand_prob_3 = np.random.rand()
               if rand_prob_3 <= mtc_m_protcount[-1]/mtc_tot_protcount[-1]:      #Remove one mutant protein
                   if mtc_m_protcount[-1] > 0:
                       mtc_m_protcount = np.append(mtc_m_protcount,np.subtract(mtc_m_protcount[-1],1))
                       mtc_tot_protcount = np.append(mtc_tot_protcount, np.add(mtc_m_protcount[-1], mtc_w_protcount))
                   t_rec.append(t)  
                   
               else:
                   if mtc_w_protcount[-1] >0:
                       mtc_w_protcount = np.append(mtc_w_protcount, np.subtract(mtc_w_protcount[-1],1))
                       mtc_tot_protcount = np.append(mtc_tot_protcount, np.add(mtc_w_protcount[-1],mtc_m_protcount))
                   t_rec.append(t)
                   
           elif prob_prod + prob_decay < rand_prob <= prob_prod + prob_decay + prob_mito: # background mitophagy events 
                  #Reset protein count array and start tracking for new mitochondria
                  mtc_wnu_protcount = np.array([]) 
                  mtc_mnu_protcount = np.array([])
                  
                           
                  #binomial distribution to split wildtype proteins
                           
                  n, p = mtc_w_protcount[-1], 1-bite_size_true #total number of proteins at split, probability of going to mtc 1 or 2
             
                  r = binom.rvs(n, p)
                           
                  mtc_wnu_protcount = np.append(mtc_wnu_protcount.astype(int),r) #add split of wildtype protein to mtc1
                  
                           
                           
                  #binomial distribution to split mutant proteins
                           
                  n, p = mtc_m_protcount[-1], 1 - bite_size_true #total number of proteins at split, probability of going to mtc 1 or 2
             
                  r = binom.rvs(n, p)
                           
                  mtc_mnu_protcount = np.append(mtc_mnu_protcount.astype(int),r) #add split of mutant protein to mtc1
                  
                           
                  mtc_totnu_protcount = np.array(np.add(mtc_wnu_protcount,mtc_mnu_protcount))
                           
                  mtc_w_protcount = mtc_wnu_protcount
                  mtc_m_protcount = mtc_mnu_protcount
                  mtc_tot_protcount = mtc_totnu_protcount
                  
                  
                  t_rec.append(t)
                  
                     
           #Update protein count in master lists    
           wcount = np.append(wcount,mtc_w_protcount[-1])
           totcount = np.append(totcount, mtc_tot_protcount[-1])
           mcount = np.append(mcount, mtc_m_protcount[-1])
                
           
           if wcount[-1] <= totcount[0]*threshold:
              t_fin.append(t)
              #wcount[-1] == 0
              break
        # if len(bite_sizes) > 0:
        #     bite_size_true_avg = sum(bite_sizes)/len(bite_sizes)
        #     bite_size_rec.append(bite_size_true_avg)
        # else:
        #     bite_size_rec.append(0)
            
        bite_size_true_avg = sum(bite_sizes)/len(bite_sizes)
        bite_size_rec.append(bite_size_true_avg)
        mito_count += 1 
        if mito_count%100==0:
            print(mito_count)
        #print(mito_count)
    #Update Master DataFrames with Master Lists
        mtcm_wcount = pd.concat([mtcm_wcount, pd.DataFrame(wcount)], axis = 1)
        mtcm_mcount = pd.concat([mtcm_mcount, pd.DataFrame(mcount)], axis = 1)
        mtcm_totcount = pd.concat([mtcm_totcount, pd.DataFrame(totcount)], axis =1)
        tm_tot = pd.concat([tm_tot, pd.DataFrame(t_rec)], axis = 1)   
    
    return mtcm_wcount,mtcm_mcount,mtcm_totcount,tm_tot, t_fin,mito_rec,decay_rec,bite_size_rec,bite_draw_rec