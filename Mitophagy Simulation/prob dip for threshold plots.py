# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:25:01 2024

@author: eshan_user
"""

import numpy as np
import matplotlib.pyplot as plt


b = np.linspace(0.05,0.5,10)
thresh = 0.25

n = np.ceil(np.log(thresh)/np.log(1-b))
n2 = np.ceil(np.log(thresh)/np.log(1-b/2))
n3 = np.ceil(np.log(thresh)/np.log(1-b/3))

plt.plot(b,n,label  = 'Mito 1',color ='r')
plt.plot(b,n2,label = 'Mito 2',color ='b')
plt.plot(b, n3, label = 'Mito 3',color ='g')
plt.legend()
plt.xlabel('Bite Size')
plt.ylabel('Number of Bites')
#plt.title('Bite counts to reach threshold')
plt.show()


ratio_12 = n/n2

ratio_13 = n/n3

plt.plot(b,ratio_12,label = 'Mito 1:Mito2',color = 'b')
plt.plot(b,ratio_13, label = 'Mito 1:Mito 3',color ='g')
plt.xlabel('Bite Size')
plt.ylabel('Ratio of bite counts')
plt.legend()
plt.show()