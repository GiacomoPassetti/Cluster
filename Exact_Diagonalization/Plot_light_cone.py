# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 11:02:25 2021

@author: giaco
"""
import matplotlib.pyplot as plt
import numpy as np

n_i_t=np.load('C:/users/giaco/Desktop/Cluster/Exact_Diagonalization/Time_occup_Exact_Omega_10J_1 g_0 Nmax_6 L_8.npy')
nit=[]

for t in range(len(n_i_t)):
    ni=[]
    for i in range(len(n_i_t[t])):
        ni.append(np.real(n_i_t[t][i]))
    nit.append(ni)
        
    
print(n_i_t[0])
plt.figure()
plt.imshow(nit[::-1],
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(0, 8, 0, 3))
plt.xlabel('site i')
plt.ylabel('time g='+str(0))

plt.colorbar().set_label('Occupancy $N$')