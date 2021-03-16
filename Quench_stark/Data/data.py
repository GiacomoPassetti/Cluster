# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 09:02:20 2021

@author: giaco
"""

import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quench_stark/Data')
import numpy as np

import matplotlib.pyplot as plt

Nmax=15
L=12
g_0=0.0
g=0
Omega  = 10
J=1
h=0.1
V=0

ID='C:/Users/giaco/Desktop/Cluster/Quench_stark/Data_12/Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g_0)



#Generate the GS from the initial Ansatz




a=[np.load(ID+'phot.npy')]
b=[np.load(ID+'occup.npy')]

gg=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,  1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]
ind=np.arange(1,13,1)
print(gg)
for g in gg:
    print('finding ground for:', g)
    ID='C:/Users/giaco/Desktop/Cluster/Quench_stark/Data_12/Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g'+str(g)
    a.append(np.load(ID+'phot.npy'))
    b.append(np.load(ID+'occup.npy'))
    

gg=[0]+list(np.arange(0.2,2,0.1))
plt.plot(gg, a, 'r+')
plt.xlabel(r'$g$')
plt.ylabel(r'$<a^{\dagger}a>$')
plt.legend(['L=12'])


"""
plt.plot(ind, b[9])
plt.xlabel(r'$g$')
plt.ylabel(r'$<c^{\dagger}c>$')
plt.legend(['L=24'])
"""