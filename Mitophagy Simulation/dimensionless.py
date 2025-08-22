# -*- coding: utf-8 -*-
"""
Created on Wed May  1 20:40:09 2024

@author: eshan_user
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def steady(km_back,f):
    kd = 1/19 - km_back*f
    return kd

def mito_func(kp,kd_p,f,t,N0):
    Nt = N0*np.exp(-kd_p*t) + (kp/(kd_p)) * (1 - np.exp(-kd_p*t))
    km_active = (kd_p - 1/19)/f
    return Nt, km_active

def pick_kd_p(kp,bite_size,half_life,N0,drop):
    kdp_try = 0.1
    kdp_int = 0.01
    while kdp_int > 0.0001:
        Nt,km_active = mito_func(kp,kdp_try,bite_size,half_life,N0)
        if Nt > N0*drop:
            kdp_try = kdp_try + kdp_int
        else:
            kdp_try = kdp_try - kdp_int
            kdp_int = kdp_int/10
            kdp_try = kdp_try + kdp_int
    return kdp_try

def eta_func(N_scale,t_scale,N0,tau,kd_p,kp):
    return  (N0/N_scale)*np.exp(-kd_p*tau*t_scale) - (kp/(kd_p*N_scale))*np.exp(-kd_p*tau*t_scale) + (kp/(kd_p*N_scale))


def eta_base(eta0,tau):
    return np.exp(-tau)*(eta0-1) + 1

eta0 = 100
taus = np.linspace(0,20,1000)
kp = (eta0/19)
bite_size = 0.25
half_life = 2
drop = 0.5
km_back = 0.01
base = 1



kd = steady(km_back,bite_size)
kp_scale = 0.5
kp_active = kp*kp_scale

etas_base = eta_base(eta0,taus)

drop1 = 0.5
kdp1 = pick_kd_p(kp_active,bite_size,half_life,eta0,drop1)
N_init1 = eta0 * (kp_active/kdp1)

N_scale1 = kp_active/kdp1
t_scale1 = 1/kdp1

#etas1 = eta_func(base,base,N_init1,taus,kdp1,kp_active)
etas1 = eta_func(N_scale1,t_scale1,N_init1,taus,kdp1,kp_active)

drop2 = 0.25
kdp2 = pick_kd_p(kp_active,bite_size,half_life,eta0,drop2)
N_init2 = eta0 * (kp_active/kdp2)

N_scale2 = kp_active/kdp2
t_scale2 = 1/kdp2

#etas2 = eta_func(base,base,N_init2,taus,kdp2,kp_active)
etas2 = eta_func(N_scale2,t_scale2,N_init2,taus,kdp2,kp_active)

drop3 = 0.75
kdp3 = pick_kd_p(kp_active,bite_size,half_life,eta0,drop3)
N_init3 = eta0 * (kp_active/kdp3)
N_scale3 = kp_active/kdp3
t_scale3 = 1/kdp3

#etas3 = eta_func(base,base,N_init3,taus,kdp3,kp_active)
etas3 = eta_func(N_scale3,t_scale3,N_init3,taus,kdp3,kp_active)


drop4 = 0.05
kdp4 = pick_kd_p(kp_active,bite_size,half_life,eta0,drop4)
N_init4 = eta0 * (kp_active/kdp4)
N_scale4 = kp_active/kdp4
t_scale4 = 1/kdp4

#etas4 = eta_func(base,base, N_init4, taus,kdp4, kp_active) 
etas4 = eta_func(N_scale4,t_scale4, N_init4, taus,kdp4, kp_active) 

fig,ax = plt.subplots(figsize = (10,10))
ax.plot(taus,etas_base,alpha = 0.3,label = 'None',lw = 15,zorder=1)
ax.plot(taus,etas4,label = r'$N_0 = $' + str(round(N_scale4,2)) + r' $t_0 = $' +str(round(t_scale4,2)),lw =5,color ='orange')
ax.plot(taus,etas2,label  = r'$N_0 = $' + str(round(N_scale2,2)) + r' $t_0 = $' +str(round(t_scale2,2)),lw = 5,color = 'g')
ax.plot(taus,etas1,label = r'$N_0 = $' + str(round(N_scale1,2)) + r' $t_0 = $' +str(round(t_scale1,2)),lw = 5,color = 'c')
ax.plot(taus,etas3,label = r'$N_0 = $' + str(round(N_scale3,2)) + r' $t_0 = $' +str(round(t_scale3,2)),lw =5,color = 'm')
ax.set_xlabel(r'$\theta$',fontsize = 26)
ax.set_ylabel(r'$\nu$',fontsize = 26)
#ax.set_title('Dimensionless protein loss varying kd\' with drop factor. kp_active = kp*' + str(kp_scale),fontsize = 30)
ax.legend(title = 'Scaling factor',fontsize = 26,title_fontsize=30)
ax.tick_params(axis='both',which='major',labelsize=24)
#ax.set_ylim(0,N0)
plt.show()


drops = np.linspace(0.1,0.75,10)
#drop0 = 0.5
etas = np.linspace(0,10000,50)
N_inits = []
N_scales = []
t_scales = []
for i in range(len(drops)):
    kdp = pick_kd_p(kp_active,bite_size,half_life,eta0,drops[i])
    print(kdp)
    N = eta0 * (kp_active/kdp)
    N_inits.append(N)
    
    N_scale = kp_active/kdp
    N_scales.append(N_scale)
    
    t_scale = 1/kdp
    t_scales.append(t_scale)


# plt.plot(drops,N_inits)
# plt.xlabel('Fraction of proteins remaining')
# plt.ylabel('Initial protein count')
# #plt.title('Required initial protein count for non-dimensionalization')
# plt.show()

plt.plot(drops,N_scales,label = r'$\nu$-scaling')
plt.plot(drops,t_scales, label = r'$\theta$-scaling')
plt.xlabel('Fraction of protein remaining')
plt.ylabel('Scaling value')
#plt.title('Scaling factors required for non-dimensionalization')
plt.legend()
plt.show()

N0 = 100
kp = (N0/19)
bite_size = 0.25
half_life = 2
drop = 0.5
km_back = 0.01
base = 1
kd = steady(km_back,bite_size)
kp_scale = 0.5
kp_active = kp*kp_scale


t = np.linspace(0,20,1000)

drop1 = 0.5
kdp1 = pick_kd_p(kp_active,bite_size,half_life,N0,drop1)
Nt1 = mito_func(kp_active,kdp1,bite_size,t,N0)[0]

drop2 = 0.25
kdp2 = pick_kd_p(kp_active,bite_size,half_life,N0,drop2)
Nt2 = mito_func(kp_active,kdp2,bite_size,t,N0)[0]

drop3 = 0.75
kdp3 = pick_kd_p(kp_active,bite_size,half_life,N0,drop3)
Nt3 = mito_func(kp_active,kdp3,bite_size,t,N0)[0]


drop4 = 0.05
kdp4 = pick_kd_p(kp_active,bite_size,half_life,N0,drop4)
Nt4 = mito_func(kp_active,kdp4,bite_size,t,N0)[0]



fig,ax = plt.subplots(figsize = (12,12))
#ax.plot(taus,etas_base,label = 'Dimensionless',lw = 6)
ax.plot(t,Nt4,label = '0.05',lw =3 ,color='orange')
ax.plot(t,Nt2,label  = '0.25',lw = 3,color='g')
ax.plot(t,Nt1,label = '0.5',lw = 3,color='c')
ax.plot(t,Nt3,label = '0.75',lw =3,color = 'm')

ax.set_xlabel('Time',fontsize = 26)
ax.set_ylabel('Protein count',fontsize = 26)
#ax.set_title('Dimensionless protein loss varying kd\' with drop factor. kp_active = kp*' + str(kp_scale),fontsize = 30)
ax.legend(title = 'Fraction of protein remaining',fontsize = 26,title_fontsize=30)
ax.tick_params(axis='both',which='major',labelsize=24)
#ax.set_ylim(0,N0)
plt.show()
