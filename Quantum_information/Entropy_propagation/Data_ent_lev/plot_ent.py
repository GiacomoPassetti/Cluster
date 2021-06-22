# -*- coding: utf-8 -*-
"""
Created on Tue May 25 08:16:02 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt
import sys


ts = np.arange(0,1.025, 0.025)
Fold = 'C:/Users/giaco/Desktop/Cluster/Quantum_information/Entropy_propagation/Data_ent_lev/'
def Data(L, hmax, level):

    a = 'Fidelty_hmax'+str(hmax)+'L_'+str(L)+'iterations_10Level_'+str(level)+'avg_1.npy'
    
    return a
L=100
lev = 5
dataset = [0.0, 0.5, 1.0, 2.0, 5.0]
plt.figure(dpi=800)
for i in dataset:
    data= np.load(Fold+Data(L, i, lev))
    
    plt.plot(ts, data)
    
plt.legend(['hmax=0.0', 'hmax=0.5','hmax=1.0','hmax=2.0','hmax=5.0'])
plt.title('Loschmidt  L='+str(L)+'  level = '+str(lev))
plt.ylabel(r'$<\psi(0)|\psi(t)>$')
plt.xlabel(r'$t$')
    sdvds