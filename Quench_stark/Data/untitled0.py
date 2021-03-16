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
L=24
g_0=0.0
g=0
Omega  = 10
J=1
h=0.1
V=0

ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g_0)



#Generate the GS from the initial Ansatz





a=np.load(ID+'phot.npy')
b=np.load(ID+'occup.npy')

gg=np.arange(0.2,2,0.1)
for g in gg:
    print('finding ground for:', g)
    ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g'+str(g)
    np.save(ID+'phot.npy', n_ph)
    np.save(ID+'occup.npy', n_f)