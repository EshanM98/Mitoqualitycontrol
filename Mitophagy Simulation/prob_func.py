# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 21:16:23 2023

@author: eshan_user
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.special as scspec
from scipy import integrate 



def size_lim_eq(size_prob,f_size,f_lim): #Distribute probability under min size equally throughout all sizes about min size - led to conclusion that some mito are too small for more than one piecemeal mitophagy event
    prob = size_prob.copy()
    mirror_prob = []
    mirror_size = []
    for i in range(len(prob)):
        if f_size[i] <= f_lim  :
            mirror_prob.append(prob[i])
            mirror_size.append(f_size[i])
    
    updated_prob = prob[len(mirror_prob):]
    p_rec = []
    p_under = sum(mirror_prob)
    p_over = sum(prob) - p_under
    for i in range(len(updated_prob)):
        p_new = updated_prob[i]/(p_over/(p_over+p_under))
        p_rec.append(p_new)
    min_size = len(f_size) - len(updated_prob)
    updated_size = f_size[min_size:]    
    return p_rec, updated_size

def init(x,prob,mu,sigma,f_lim): #Initializes and normalizes size distribution
    gauss = prob*np.exp((-1/2)*((x-mu)**2)/(sigma**2))/(sigma*np.sqrt(2*np.pi))
    gauss_lim,size_fix = size_lim_eq(gauss,x,f_lim)
    y = integrate.cumulative_trapezoid(gauss_lim,size_fix)
    gauss_fix = np.asarray(gauss_lim)/(y[-1])
    y_int = integrate.cumulative_trapezoid(gauss_fix,size_fix)
    return gauss_fix,size_fix,y_int


def choose_size(size,prob):
    y_int = integrate.cumulative_trapezoid(prob,size)
    rand_prob = np.random.uniform(0,y_int[-1])
    idx = np.abs(np.asarray(y_int)-rand_prob).argmin()
    size_choose = size[idx]
    return size_choose

def poisson(k,lamb):
    return lamb**k * np.exp(-lamb)/scspec.factorial(k)

def norm(k,lamb):
    return (1/np.sqrt(2*np.pi*lamb)) * np.exp(-0.5*((k-lamb)/np.sqrt(lamb)**2))

        
