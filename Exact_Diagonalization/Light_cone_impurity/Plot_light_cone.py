# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 11:02:25 2021

@author: giaco
"""
import matplotlib.pyplot as plt
import numpy as np

n_i_t=np.load('C:/users/giaco/Desktop/Cluster/Exact_Diagonalization/nit_0-5_Omega_10J_1 g_0 Nmax_6 L_8.npy')


        
    
plt.figure()
plt.imshow(n_i_t[::-1],'plasma',
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(0, 8, 0, 5))
plt.xlabel('site i')
plt.ylabel('time g='+str(1))

plt.colorbar().set_label('Occupancy $N$')