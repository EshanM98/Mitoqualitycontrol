# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 11:06:32 2023

@author: eshan_user
"""

import numpy as np

def binning(bins,binthreshold,prevgen_all,nextgen_all):
    
    
    group_x = np.linspace(0,100,bins+1)
    groupsprev_avg = []
    groupsnext_avg = []
    groups_prev = []
    groups_next = []
    for i in range(len(group_x)-1):
        bincount = 0
        group_prev = []
        group_next = []
        for j in range(len(prevgen_all)):
            if group_x[i] <= prevgen_all[j] < group_x[i+1]:
                bincount += 1
                group_prev.append(prevgen_all[j])
                group_next.append(nextgen_all[j])
        if bincount >= binthreshold:
            groupsprev_avg.append(sum(group_prev)/len(group_prev))
            groupsnext_avg.append(sum(group_next)/len(group_next))
            
        groups_prev.append(group_prev)
        groups_next.append(group_next)
    return groupsprev_avg,groupsnext_avg,groups_prev,groups_next

#Degradation model functions


def deg2A(f0): #Function takes initial mutation level
    n = np.arange(0,6,1) #number of mtDNA in a mtc
    f = f0 #Prob that mtDNA copy is mutant
    p_m = [0.739,0.226,0.035] #prob of having some number of nucleoids in mtc
    p_w = [0.64,0.36] # prob of having some number of DNA copies in nucleoid
    p_n = np.zeros_like(n).astype(np.float32) #Prob of having some number of mtDNA in mtc
    
    #Initialize p_n probs
    p_n[0] = p_m[0]*p_w[0]  #Prob of having 1 mtDNA in mtc
    p_n[1] = p_m[0]*p_w[1] + p_m[1]*p_w[0]**2 #prob of having 2 mtDNA in mtc
    p_n[2] = p_m[1]*2*p_w[1]*p_w[0] + p_m[2]*p_w[0]**3 # prob of having 3 mtDNA in mtc
    p_n[3] = p_m[1]*p_w[1]**2 + p_m[2]*3*p_w[1]*p_w[0]**2 #prob of having 4 mtDNA in mtc
    p_n[4] = p_m[2]*3*p_w[1]**2*p_w[0] #prob of having 5 mtDNA in mtc
    p_n[5] = p_m[2]*p_w[1]**3 #prob of having 6 mtDNA in mtc
    
    
    t2 = 1* p_n[0] + 2*p_n[1] + 3*p_n[2] + 4*p_n[3] + 5*p_n[4] + 6*p_n[5]
    
    
    t1 =  1* p_n[0]*f + 2*p_n[1]*f**2 + 3*p_n[2]*f**3 + 4*p_n[3]*f**4 + 5*p_n[4]*f**5 + 6*p_n[5]*f**6
    f_next = (f*t2 - t1)/(t2-t1)
        
    return f_next,t2,t1 #Returns the mutation level in the next generation

def deg1(f0): #Function takes initial mutation level
    n = np.arange(0,8,1) #number of mtDNA in a mtc
    f = f0 #Prob that mtDNA copy is mutant


    p_m = [0.404,0.36,0.174,0.062] #prob of having some number of nucleoids in mtc
    p_w = [0.64,0.36] # prob of having some number of DNA copies in nucleoid
    p_n = np.zeros_like(n).astype(np.float32) #Prob of having some number of mtDNA in mtc

    #Initialize p_n probs
    p_n[0] = p_m[0]*p_w[0]  #Prob of having 1 mtDNA in mtc
    p_n[1] = p_m[0]*p_w[1] + p_m[1]*p_w[0]**2 #prob of having 2 mtDNA in mtc
    p_n[2] = p_m[1]*2*p_w[1]*p_w[0] + p_m[2]*p_w[0]**3 # prob of having 3 mtDNA in mtc
    p_n[3] = p_m[1]*p_w[1]**2 + p_m[2]*3*p_w[1]*p_w[0]**2 + p_m[3]*p_w[0]**4#prob of having 4 mtDNA in mtc
    p_n[4] = p_m[2]*3*p_w[1]**2*p_w[0] + p_m[3]*4*p_w[1]*p_w[0]**3#prob of having 5 mtDNA in mtc
    p_n[5] = p_m[2]*p_w[1]**3 +p_m[3]*6*p_w[1]**2*p_w[0]**2#prob of having 6 mtDNA in mtc
    p_n[6] = p_m[3]*4*p_w[0]*p_w[1]**3
    p_n[7] = p_m[3]*p_w[1]**4


    t2 = 1* p_n[0] + 2*p_n[1] + 3*p_n[2] + 4*p_n[3] + 5*p_n[4] + 6*p_n[5] + 7*p_n[6] + 8*p_n[7] 
    t1 =  1* p_n[0]*f + 2*p_n[1]*f**2 + 3*p_n[2]*f**3 + 4*p_n[3]*f**4 + 5*p_n[4]*f**5 + 6*p_n[5]*f**6 + 7*p_n[6]*f**7 + 8*p_n[7]*f**8
    f_next = (f*t2 - t1)/(t2-t1)
        
    return f_next,t2,t1 #Returns mutation level of next generation

def inter_reg(f0,p):
    f_reg1 = deg1(f0)
    f_reg2 = deg2A(f0)
    
    f_inter = p*f_reg1[0] + (1-p)*f_reg2[0]
    return f_inter

#Square difference for degradation model
def sse(prevgen,nextgen,p): #Weighed by dividing SD by previous mutant gen (Sum of Squared Errors)
    sd_reg1 = []
    for i in range(len(prevgen)):   
        y0 = deg1(prevgen[i]/100)[0]
        y = nextgen[i]/100
        diff = ((y-y0)/(prevgen[i]/100))**2
        sd_reg1.append(diff)
    sse_r1 = sum(sd_reg1)
    
    sd_reg2 = []
    for i in range(len(prevgen)):
        y0 = deg2A(prevgen[i]/100)[0]
        y = nextgen[i]/100
        diff = ((y-y0)/(prevgen[i]/100))**2
        sd_reg2.append(diff)
    sse_r2 = sum(sd_reg2)
    
    sd_int = []
    for i in range(len(prevgen)):
        y0 = inter_reg(prevgen[i]/100,p)
        y = nextgen[i]/100
        diff = ((y-y0)/(prevgen[i]/100))**2
        sd_int.append(diff)
    sse_int = sum(sd_int)
    
    
    return sse_r1, sse_r2,sse_int

def deg2A_onerem(f0): #Function takes initial mutation level
    n = np.arange(0,6,1) #number of mtDNA in a mtc
    f = f0 #Prob that mtDNA copy is mutant
    p_m = [0.739,0.226,0.035] #prob of having some number of nucleoids in mtc
    p_w = [0.64,0.36] # prob of having some number of DNA copies in nucleoid
    p_n = np.zeros_like(n).astype(np.float32) #Prob of having some number of mtDNA in mtc
    
    #Initialize p_n probs
    p_n[0] = p_m[0]*p_w[0]  #Prob of having 1 mtDNA in mtc
    p_n[1] = p_m[0]*p_w[1] + p_m[1]*p_w[0]**2 #prob of having 2 mtDNA in mtc
    p_n[2] = p_m[1]*2*p_w[1]*p_w[0] + p_m[2]*p_w[0]**3 # prob of having 3 mtDNA in mtc
    p_n[3] = p_m[1]*p_w[1]**2 + p_m[2]*3*p_w[1]*p_w[0]**2 #prob of having 4 mtDNA in mtc
    p_n[4] = p_m[2]*3*p_w[1]**2*p_w[0] #prob of having 5 mtDNA in mtc
    p_n[5] = p_m[2]*3*p_w[1]**3 #prob of having 6 mtDNA in mtc
    
    
    t2 = 1* p_n[0] + 2*p_n[1] + 3*p_n[2] + 4*p_n[3] + 5*p_n[4] + 6*p_n[5] #tot number of dna copies
    
    
    t1 = 1*p_m[0]*p_w[0]*f + 2*p_m[0]*p_w[1]*f**2 #number of purely mutant dna copies
    f_next = (f*t2 - t1)/(t2-t1)
        
    return f_next #Returns the mutation level in the next generation

def deg2A_tworem(f0): #Function takes initial mutation level
    n = np.arange(0,6,1) #number of mtDNA in a mtc
    f = f0 #Prob that mtDNA copy is mutant
    p_m = [0.739,0.211,0.05] #prob of having some number of nucleoids in mtc
    p_w = [0.64,0.36] # prob of having some number of DNA copies in nucleoid
    p_n = np.zeros_like(n).astype(np.float32) #Prob of having some number of mtDNA in mtc
    
    #Initialize p_n probs
    p_n[0] = p_m[0]*p_w[0]  #Prob of having 1 mtDNA in mtc
    p_n[1] = p_m[0]*p_w[1] + p_m[1]*p_w[0]**2 #prob of having 2 mtDNA in mtc
    p_n[2] = p_m[1]*2*p_w[1]*p_w[0] + p_m[2]*p_w[0]**3 # prob of having 3 mtDNA in mtc
    p_n[3] = p_m[1]*p_w[1]**2 + p_m[2]*3*p_w[1]*p_w[0]**2 #prob of having 4 mtDNA in mtc
    p_n[4] = p_m[2]*3*p_w[1]**2*p_w[0] #prob of having 5 mtDNA in mtc
    p_n[5] = p_m[2]*3*p_w[1]**3 #prob of having 6 mtDNA in mtc
    
    
    t2 = 1* p_n[0] + 2*p_n[1] + 3*p_n[2] + 4*p_n[3] + 5*p_n[4] + 6*p_n[5]
    
    t1 = 1*p_m[0]*p_w[0]*f + 2*p_m[0]*p_w[1]*f**2 + 2*p_m[1]*p_w[0]**2*f**2 + 3*2*p_m[1]*p_w[1]*p_w[0]*f**3 + 4*p_m[1]*p_w[1]**2*f**4
    f_next = (f*t2 - t1)/(t2-t1)
        
    return f_next #Returns the mutation level in the next generation


def deg2A_int(f0,p):
    f_onerem = deg2A_onerem(f0)
    f_tworem = deg2A_tworem(f0)
    
    f_inter = p*f_onerem + (1-p)*f_tworem
    return f_inter

def deg1_onerem(f0): #Function takes initial mutation level
    n = np.arange(0,8,1) #number of mtDNA in a mtc
    f = f0 #Prob that mtDNA copy is mutant


    p_m = [0.404,0.36,0.174,0.062] #prob of having some number of nucleoids in mtc
    p_w = [0.64,0.36] # prob of having some number of DNA copies in nucleoid
    p_n = np.zeros_like(n).astype(np.float32) #Prob of having some number of mtDNA in mtc

    #Initialize p_n probs
    p_n[0] = p_m[0]*p_w[0]  #Prob of having 1 mtDNA in mtc
    p_n[1] = p_m[0]*p_w[1] + p_m[1]*p_w[0]**2 #prob of having 2 mtDNA in mtc
    p_n[2] = p_m[1]*2*p_w[1]*p_w[0] + p_m[2]*p_w[0]**3 # prob of having 3 mtDNA in mtc
    p_n[3] = p_m[1]*p_w[1]**2 + p_m[2]*3*p_w[1]*p_w[0]**2 + p_m[3]*p_w[0]**4#prob of having 4 mtDNA in mtc
    p_n[4] = p_m[2]*3*p_w[1]**2*p_w[0] + p_m[3]*4*p_w[1]*p_w[0]**3#prob of having 5 mtDNA in mtc
    p_n[5] = p_m[2]*3*p_w[1]**3 +p_m[3]*6*p_w[1]**2*p_w[0]**2#prob of having 6 mtDNA in mtc
    p_n[6] = p_m[3]*4*p_w[0]*p_w[1]**3
    p_n[7] = p_m[3]*p_w[1]**4


    t2 = 1* p_n[0] + 2*p_n[1] + 3*p_n[2] + 4*p_n[3] + 5*p_n[4] + 6*p_n[5] + 7*p_n[6] + 8*p_n[7] 
    t1 = 1*p_m[0]*p_w[0]*f + 2*p_m[0]*p_w[1]*f**2
    f_next = (f*t2 - t1)/(t2-t1)
        
    return f_next #Returns mutation level of next generation

def deg1_tworem(f0): #Function takes initial mutation level
    n = np.arange(0,8,1) #number of mtDNA in a mtc
    f = f0 #Prob that mtDNA copy is mutant


    p_m = [0.404,0.36,0.174,0.062] #prob of having some number of nucleoids in mtc
    p_w = [0.64,0.36] # prob of having some number of DNA copies in nucleoid
    p_n = np.zeros_like(n).astype(np.float32) #Prob of having some number of mtDNA in mtc

    #Initialize p_n probs
    p_n[0] = p_m[0]*p_w[0]  #Prob of having 1 mtDNA in mtc
    p_n[1] = p_m[0]*p_w[1] + p_m[1]*p_w[0]**2 #prob of having 2 mtDNA in mtc
    p_n[2] = p_m[1]*2*p_w[1]*p_w[0] + p_m[2]*p_w[0]**3 # prob of having 3 mtDNA in mtc
    p_n[3] = p_m[1]*p_w[1]**2 + p_m[2]*3*p_w[1]*p_w[0]**2 + p_m[3]*p_w[0]**4#prob of having 4 mtDNA in mtc
    p_n[4] = p_m[2]*3*p_w[1]**2*p_w[0] + p_m[3]*4*p_w[1]*p_w[0]**3#prob of having 5 mtDNA in mtc
    p_n[5] = p_m[2]*3*p_w[1]**3 +p_m[3]*6*p_w[1]**2*p_w[0]**2#prob of having 6 mtDNA in mtc
    p_n[6] = p_m[3]*4*p_w[0]*p_w[1]**3
    p_n[7] = p_m[3]*p_w[1]**4


    t2 = 1* p_n[0] + 2*p_n[1] + 3*p_n[2] + 4*p_n[3] + 5*p_n[4] + 6*p_n[5] + 7*p_n[6] + 8*p_n[7] 
    t1 = 1*p_m[0]*p_w[0]*f + 2*p_m[0]*p_w[1]*f**2 + 2*p_m[1]*p_w[0]**2*f**2 + 3*2*p_m[1]*p_w[1]*p_w[0]*f**3 + 4*p_m[1]*p_w[1]**2*f**4
    f_next = (f*t2 - t1)/(t2-t1)
        
    return f_next #Returns mutation level of next generation

def deg1_threerem(f0): #Function takes initial mutation level
    n = np.arange(0,8,1) #number of mtDNA in a mtc
    f = f0 #Prob that mtDNA copy is mutant


    p_m = [0.404,0.36,0.174,0.062] #prob of having some number of nucleoids in mtc
    p_w = [0.64,0.36] # prob of having some number of DNA copies in nucleoid
    p_n = np.zeros_like(n).astype(np.float32) #Prob of having some number of mtDNA in mtc

    #Initialize p_n probs
    p_n[0] = p_m[0]*p_w[0]  #Prob of having 1 mtDNA in mtc
    p_n[1] = p_m[0]*p_w[1] + p_m[1]*p_w[0]**2 #prob of having 2 mtDNA in mtc
    p_n[2] = p_m[1]*2*p_w[1]*p_w[0] + p_m[2]*p_w[0]**3 # prob of having 3 mtDNA in mtc
    p_n[3] = p_m[1]*p_w[1]**2 + p_m[2]*3*p_w[1]*p_w[0]**2 + p_m[3]*p_w[0]**4#prob of having 4 mtDNA in mtc
    p_n[4] = p_m[2]*3*p_w[1]**2*p_w[0] + p_m[3]*4*p_w[1]*p_w[0]**3#prob of having 5 mtDNA in mtc
    p_n[5] = p_m[2]*3*p_w[1]**3 +p_m[3]*6*p_w[1]**2*p_w[0]**2#prob of having 6 mtDNA in mtc
    p_n[6] = p_m[3]*4*p_w[0]*p_w[1]**3
    p_n[7] = p_m[3]*p_w[1]**4


    t2 = 1* p_n[0] + 2*p_n[1] + 3*p_n[2] + 4*p_n[3] + 5*p_n[4] + 6*p_n[5] + 7*p_n[6] + 8*p_n[7] 
    t1 = 1*p_m[0]*p_w[0]*f + 2*p_m[0]*p_w[1]*f**2 + 2*p_m[1]*p_w[0]**2*f**2 + 3*2*p_m[1]*p_w[1]*p_w[0]*f**3 + 3*p_m[2]*p_w[0]**3*f**3 + 4*p_m[1]*p_w[1]**2*f**4 + 4*p_m[2]*3*p_w[1]*p_w[0]**2*f**4 + 5*p_m[2]*3*p_w[1]**2*p_w[0]*f**5 + 6*p_m[2]*3*p_w[1]**3*f**6
    f_next = (f*t2 - t1)/(t2-t1)
        
    return f_next #Returns mutation level of next generation


#Replication model

def rep_comp(f1,f2,n,lamb):
    return f1*140*lamb**n + f2*140*2**n + (1-f1-f2)*140*2**n

def rep_noncomp(f,n,lamb):
    return f*140*lamb**n + (1-f)*140*2**n

def rep_nocomp(f0,lamb,mtc_fin):
    f = f0
    
    #Determine value of n for a given lambda value
    ntry = 1
    n_int = 0.1
    while n_int > 0.0001:
        mtc_count = rep_noncomp(f,ntry,lamb) 
            
        if mtc_count < mtc_fin:
            ntry = ntry + n_int
        else:
            ntry = ntry - n_int
            n_int = n_int/10
            ntry = ntry + n_int
        
    n = ntry
        
    f_nextgen_nocomp = (f*lamb**n)/(f*lamb**n + (1-f)*2**n)
    return f_nextgen_nocomp

def rep_r1(f0,lamb,mtc_fin):
    p_m = [0.404,0.36,0.174,0.062] #prob of having some number of nucleoids in mtc
    p_w = [0.64,0.36] # prob of having some number of DNA copies in nucleoid
    p_n = np.zeros(8).astype(np.float32) #Prob of having some number of mtDNA in mtc
    
    
    def pure_mutr1(f,p_n):
        t1 = p_n[0] * f + 2*p_n[1] * f**2 + 3*p_n[2] * f**3 + 4*p_n[3] * f**4 + 5*p_n[4] * f**5 + 6*p_n[5] * f**6 + 7*p_n[6] * f**7 + 8*p_n[7] * f**8
        return t1
    
    def f_r1(p_m,p_w,p_n):
        #Initialize p_n probs
        p_n[0] = p_m[0]*p_w[0]  #Prob of having 1 mtDNA in mtc
        p_n[1] = p_m[0]*p_w[1] + p_m[1]*p_w[0]**2 #prob of having 2 mtDNA in mtc
        p_n[2] = p_m[1]*2*p_w[1]*p_w[0] + p_m[2]*p_w[0]**3 # prob of having 3 mtDNA in mtc
        p_n[3] = p_m[1]*p_w[1]**2 + p_m[2]*3*p_w[1]*p_w[0]**2 + p_m[3]*p_w[0]**4#prob of having 4 mtDNA in mtc
        p_n[4] = p_m[2]*3*p_w[1]**2*p_w[0] + p_m[3]*4*p_w[1]*p_w[0]**3#prob of having 5 mtDNA in mtc
        p_n[5] = p_m[2]*p_w[1]**3 +p_m[3]*6*p_w[1]**2*p_w[0]**2#prob of having 6 mtDNA in mtc
        p_n[6] = p_m[3]*4*p_w[0]*p_w[1]**3
        p_n[7] = p_m[3]*p_w[1]**4
        
        
        t2 = 1* p_n[0] + 2*p_n[1] + 3*p_n[2] + 4*p_n[3] + 5*p_n[4] + 6*p_n[5] + 7*p_n[6] + 8*p_n[7] 
        return p_n,t2
    
    f = f0
    
    f_iso = pure_mutr1(f,f_r1(p_m,p_w,p_n)[0])/f_r1(p_m,p_w,p_n)[1]
    f_noniso = f - f_iso
    ntry = 1
    n_int = 0.1
    while n_int > 0.0001:
        mtc_count = rep_comp(f_iso,f_noniso,ntry,lamb)
            
        if mtc_count < mtc_fin:
            ntry = ntry + n_int
        else:
            ntry = ntry - n_int
            n_int = n_int/10
            ntry = ntry + n_int
    n = ntry
    f_nextgen_comp = (f_iso*lamb**n + f_noniso*2**n)/(f_iso*lamb**n+f_noniso*2**n+(1-f_noniso-f_iso)*2**n)
    
    return f_nextgen_comp,n
def rep_r2(f0,lamb,mtc_fin):
    
    p_m = [0.739,0.226,0.035] #prob of having some number of nucleoids in mtc
    p_w = [0.64,0.36] # prob of having some number of DNA copies in nucleoid
    p_n = np.zeros(6).astype(np.float32) #Prob of having some number of mtDNA in mtc
    
    def pure_mutr2(f,p_n):
        t1 = p_n[0] * f + 2*p_n[1] * f**2 + 3*p_n[2] * f**3 + 4*p_n[3] * f**4 + 5*p_n[4] * f**5 + 6*p_n[5] * f**6
        return t1
    
    def f_r2(p_m,p_w,p_n):
        #Initialize p_n probs
        p_n[0] = p_m[0]*p_w[0]  #Prob of having 1 mtDNA in mtc
        p_n[1] = p_m[0]*p_w[1] + p_m[1]*p_w[0]**2 #prob of having 2 mtDNA in mtc
        p_n[2] = p_m[1]*2*p_w[1]*p_w[0] + p_m[2]*p_w[0]**3 # prob of having 3 mtDNA in mtc
        p_n[3] = p_m[1]*p_w[1]**2 + p_m[2]*3*p_w[1]*p_w[0]**2 #prob of having 4 mtDNA in mtc
        p_n[4] = p_m[2]*3*p_w[1]**2*p_w[0] #prob of having 5 mtDNA in mtc
        p_n[5] = p_m[2]*3*p_w[1]**3 #prob of having 6 mtDNA in mtc
        
        
        t2 = 1* p_n[0] + 2*p_n[1] + 3*p_n[2] + 4*p_n[3] + 5*p_n[4] + 6*p_n[5] #total number of mtDNA copies
        return p_n, t2
    
    f = f0
        
    f_iso = pure_mutr2(f,f_r2(p_m,p_w,p_n)[0])/f_r2(p_m,p_w,p_n)[1]
    f_noniso = f - f_iso
    ntry = 1
    n_int = 0.1
    while n_int > 0.0001:
        mtc_count = rep_comp(f_iso,f_noniso,ntry,lamb)
            
        if mtc_count < mtc_fin:
            ntry = ntry + n_int
        else:
            ntry = ntry - n_int
            n_int = n_int/10
            ntry = ntry + n_int
    n = ntry
    f_nextgen_comp = (f_iso*lamb**n + f_noniso*2**n)/(f_iso*lamb**n+f_noniso*2**n+(1-f_noniso-f_iso)*2**n)
    return f_nextgen_comp,n

#Square differences for replication model
def sse_rep_nocomp(lamb,prevgen,nextgen,mtc_fin):
    square_diff_noncomp = []
    for i in range(len(prevgen)):
        y0 = rep_nocomp(prevgen[i]/100,lamb,mtc_fin)
        y = nextgen[i]/100
        diff = ((y-y0)/(prevgen[i]/100))**2
        square_diff_noncomp.append(diff)
    sse_noncomp = sum(square_diff_noncomp)
    return sse_noncomp

def sse_rep_r1(lamb,prevgen,nextgen,mtc_fin):
    square_diff_comp = []
    for i in range(len(prevgen)):
        y0 = rep_r1(prevgen[i]/100,lamb,mtc_fin)
        y = nextgen[i]/100
        diff = ((y-y0)/(prevgen[i]/100))**2
        square_diff_comp.append(diff)
    sse_comp = sum(square_diff_comp)
    return sse_comp

def sse_rep_r2(lamb,prevgen,nextgen,mtc_fin):
    square_diff_comp = []
    for i in range(len(prevgen)):
        y0 = rep_r2(prevgen[i]/100,lamb,mtc_fin)
        y = nextgen[i]/100
        diff = ((y-y0)/(prevgen[i]/100))**2
        square_diff_comp.append(diff)
    sse_comp = sum(square_diff_comp)
    return sse_comp



#Alternative rep models (based on number of replications)
def rep_r1_alt(f0,lamb,n):
    p_m = [0.404,0.36,0.174,0.062] #prob of having some number of nucleoids in mtc
    p_w = [0.64,0.36] # prob of having some number of DNA copies in nucleoid
    p_n = np.zeros(8).astype(np.float32) #Prob of having some number of mtDNA in mtc
    
    
    def pure_mutr1(f,p_n):
        t1 = p_n[0] * f + 2*p_n[1] * f**2 + 3*p_n[2] * f**3 + 4*p_n[3] * f**4 + 5*p_n[4] * f**5 + 6*p_n[5] * f**6 + 7*p_n[6] * f**7 + 8*p_n[7] * f**8
        return t1
    
    def f_r1(p_m,p_w,p_n):
        #Initialize p_n probs
        p_n[0] = p_m[0]*p_w[0]  #Prob of having 1 mtDNA in mtc
        p_n[1] = p_m[0]*p_w[1] + p_m[1]*p_w[0]**2 #prob of having 2 mtDNA in mtc
        p_n[2] = p_m[1]*2*p_w[1]*p_w[0] + p_m[2]*p_w[0]**3 # prob of having 3 mtDNA in mtc
        p_n[3] = p_m[1]*p_w[1]**2 + p_m[2]*3*p_w[1]*p_w[0]**2 + p_m[3]*p_w[0]**4#prob of having 4 mtDNA in mtc
        p_n[4] = p_m[2]*3*p_w[1]**2*p_w[0] + p_m[3]*4*p_w[1]*p_w[0]**3#prob of having 5 mtDNA in mtc
        p_n[5] = p_m[2]*p_w[1]**3 +p_m[3]*6*p_w[1]**2*p_w[0]**2#prob of having 6 mtDNA in mtc
        p_n[6] = p_m[3]*4*p_w[0]*p_w[1]**3
        p_n[7] = p_m[3]*p_w[1]**4
        
        
        t2 = 1* p_n[0] + 2*p_n[1] + 3*p_n[2] + 4*p_n[3] + 5*p_n[4] + 6*p_n[5] + 7*p_n[6] + 8*p_n[7] 
        return p_n,t2
    
    f = f0
    
        
    f_nextgen_nocomp = (f*lamb**n)/(f*lamb**n + (1-f)*2**n)
    
    fin_nocomp = rep_noncomp(f,n,lamb)
    f_iso = pure_mutr1(f,f_r1(p_m,p_w,p_n)[0])/f_r1(p_m,p_w,p_n)[1]
    
    f_noniso = f - f_iso
    
    f_nextgen_comp = (f_iso*lamb**n + f_noniso*2**n)/(f_iso*lamb**n+f_noniso*2**n+(1-f_noniso-f_iso)*2**n)
    fin_compr1 = rep_comp(f_iso,f_noniso,n,lamb)
    return f_nextgen_nocomp, f_nextgen_comp,fin_nocomp,fin_compr1

def rep_r2_alt(f0,lamb,n):
    
    p_m = [0.739,0.226,0.035] #prob of having some number of nucleoids in mtc
    p_w = [0.64,0.36] # prob of having some number of DNA copies in nucleoid
    p_n = np.zeros(6).astype(np.float32) #Prob of having some number of mtDNA in mtc
    
    def pure_mutr2(f,p_n):
        t1 = p_n[0] * f + 2*p_n[1] * f**2 + 3*p_n[2] * f**3 + 4*p_n[3] * f**4 + 5*p_n[4] * f**5 + 6*p_n[5] * f**6
        return t1
    
    def f_r2(p_m,p_w,p_n):
        #Initialize p_n probs
        p_n[0] = p_m[0]*p_w[0]  #Prob of having 1 mtDNA in mtc
        p_n[1] = p_m[0]*p_w[1] + p_m[1]*p_w[0]**2 #prob of having 2 mtDNA in mtc
        p_n[2] = p_m[1]*2*p_w[1]*p_w[0] + p_m[2]*p_w[0]**3 # prob of having 3 mtDNA in mtc
        p_n[3] = p_m[1]*p_w[1]**2 + p_m[2]*3*p_w[1]*p_w[0]**2 #prob of having 4 mtDNA in mtc
        p_n[4] = p_m[2]*3*p_w[1]**2*p_w[0] #prob of having 5 mtDNA in mtc
        p_n[5] = p_m[2]*3*p_w[1]**3 #prob of having 6 mtDNA in mtc
        
        
        t2 = 1* p_n[0] + 2*p_n[1] + 3*p_n[2] + 4*p_n[3] + 5*p_n[4] + 6*p_n[5] #total number of mtDNA copies
        return p_n, t2
    
    f = f0
        
        
    f_nextgen_nocomp = (f*lamb**n)/(f*lamb**n + (1-f)*2**n)
    fin_noncomp = rep_noncomp(f, n, lamb)
    
    f_iso = pure_mutr2(f,f_r2(p_m,p_w,p_n)[0])/f_r2(p_m,p_w,p_n)[1]
    f_noniso = f - f_iso
    
    f_nextgen_comp = (f_iso*lamb**n + f_noniso*2**n)/(f_iso*lamb**n+f_noniso*2**n+(1-f_noniso-f_iso)*2**n)
    fin_compr2 = rep_comp(f_iso, f_noniso, n, lamb)
    return f_nextgen_nocomp, f_nextgen_comp,fin_noncomp,fin_compr2

def sse_rep_r1_alt(lamb,prevgen,nextgen,n):
    square_diff_noncomp = []
    for i in range(len(prevgen)):
        y0 = rep_r1_alt(prevgen[i]/100,lamb,n)[0]
        y = nextgen[i]/100
        diff = ((y0-y)/(prevgen[i]/100))**2
        square_diff_noncomp.append(diff)
    sse_noncomp = sum(square_diff_noncomp)
    
    square_diff_comp = []
    for i in range(len(prevgen)):
        y0 = rep_r1_alt(prevgen[i]/100,lamb,n)[1]
        y = nextgen[i]/100
        diff = ((y0-y)/(prevgen[i]/100))**2
        square_diff_comp.append(diff)
    sse_comp = sum(square_diff_comp)
    return sse_noncomp, sse_comp

def sse_rep_r2_alt(lamb,prevgen,nextgen,n):
    square_diff_noncomp = []
    for i in range(len(prevgen)):
        y0 = rep_r2_alt(prevgen[i]/100,lamb,n)[0]
        y = nextgen[i]/100
        diff = ((y0-y)/(prevgen[i]/100))**2
        square_diff_noncomp.append(diff)
    sse_noncomp = sum(square_diff_noncomp)
    
    square_diff_comp = []
    for i in range(len(prevgen)):
        y0 = rep_r2_alt(prevgen[i]/100,lamb,n)[1]
        y = nextgen[i]/100
        diff = ((y0-y)/(prevgen[i]/100))**2
        square_diff_comp.append(diff)
    sse_comp = sum(square_diff_comp)
    return sse_noncomp, sse_comp

