# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 12:23:05 2024

@author: eshan_user
"""

import numpy as np
import matplotlib.pyplot as plt
from prob_func import *

## Initialize size distribution
# a1 = height of peak
#mu =  center of peak

size = 10000
nuc_prob = [0.739,0.211,0.05]
nuc_count = [1,2,3]
x= np.linspace(0.01,5,size)
f_lim = 0.75
sigma = [0.25]
for i in range(len(sigma)):

    
    gauss1_norm,size1,y_norm1 = init(x,nuc_prob[0],nuc_count[0],sigma[i],f_lim)
    gauss2_norm,size2,y_norm2 = init(x,nuc_prob[1],nuc_count[1],sigma[i],f_lim)
    gauss3_norm,size3,y_norm3 = init(x,nuc_prob[2],nuc_count[2],sigma[i],f_lim)
    
   
    fig,ax = plt.subplots()
    
    
    ax.plot(size1,gauss1_norm,label = '1',color = 'r')
    ax.plot(size2,gauss2_norm, label = '2',color = 'b')
    ax.plot(size3,gauss3_norm,label = '3',color = 'g')
    
    
    ax.bar(nuc_count,nuc_prob, alpha = 0.5)
    #ax.hist(nuc_prob,bins = 3)
    ax.set_xlabel('Nucleoid count')
    ax.set_xticks(np.arange(0,6,step = 1))
    ax.set_ylabel('Probability Density')
    #ax.set_title('Nucleoid count distribution')
    ax.legend(title = 'Nucleoid count')
    plt.show()




