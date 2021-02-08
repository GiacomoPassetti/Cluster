# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 11:02:25 2021

@author: giaco
"""
import matplotlib.pyplot as plt
import numpy as np



n_i_t=np.load('D:/Data_Spinless_boson/Coherent_state/LC_coherent_L120_g1.0_Omega_5.0displacement_1.0dt_0.01nit.npy')

print(len(n_i_t))
        
    
plt.figure()
plt.imshow(n_i_t[::-1],'plasma',
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(0, 120, 0, 3))
plt.xlabel('site i')
plt.ylabel('time t')

plt.colorbar().set_label('Occupancy $N$')


"""
sv=np.load('D:/Data_Spinless_boson/Coherent_state/LC_coherent_L120_g1.0_Omega_5.0displacement_1.0dt_0.01X(t).npy')

t=np.arange(0,3, 0.01)
plt.plot(t,sv)
plt.xlabel('t')
plt.ylabel('SV discarded')
"""

"""
xt=np.load('D:/Data_Spinless_boson/Coherent_state/LC_coherent_L120_g2.0_Omega_5.0displacement_0.25dt_0.01X(t).npy')
t=np.arange(0,3, 0.01)
plt.plot(t,xt)
plt.xlabel('t')
plt.ylabel('<x(t)>')
"""