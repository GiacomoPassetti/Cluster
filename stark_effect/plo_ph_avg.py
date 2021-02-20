# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 10:20:47 2021

@author: giaco
"""
import numpy as np
import matplotlib.pyplot as plt 



gsn=list(np.arange(0, 4, 0.2))
gs=list(np.arange(0,4,0.2))
n1=np.load('C:/Users/giaco/Desktop/Cluster/Stark_effect/Average_boson_half_filling_starkOmega_0.5J_1 g_0 Nmax_8 L_6.npy')
n2=np.load('C:/Users/giaco/Desktop/Cluster/Stark_effect/Average_boson_half_filling_starkOmega_0.5J_1 g_0 Nmax_8 L_8.npy')


plt.plot(gs, n1, 'bo', gs, n2, 'r+')
plt.ylabel(r'$<N>$')
plt.xlabel('g')
plt.legend([r'$L=6$', r'$L=8$'])
plt.show

