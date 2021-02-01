# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 18:44:56 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


"""
TEBD=np.load('C:/users/giaco/Desktop/Cluster/Exact_Diagonalization/Data/N_avg_bosEXACT_TEBDPsi_GS_Nmax_4L_8g_1.8Omega_1J_1h_0V_0.npy')
Exact=np.load('C:/users/giaco/Desktop/Cluster/Exact_Diagonalization/Data/N_avg_bosEXACTOmega_1J_1 g_1 Nmax_4 L_8.npy')
gtebd=np.arange(0,2,0.2)
gexact=np.arange(0,2,0.05)

plt.plot(gtebd, TEBD, 'b', gexact, Exact, 'r--')
plt.xlabel(r'$g$')
plt.ylabel(r'$<N_avg>_{GS}$')
red_patch = mpatches.Patch(color='red', label='Exact Diag')

blue_patch = mpatches.Patch(color='blue', label='Imaginary TE')
plt.legend(handles=[blue_patch, red_patch] )
plt.show
"""

stark=np.load('C:/users/giaco/Desktop/Cluster/Exact_Diagonalization/Data/N_avg_bosEXACT_StarkOmega_1J_1 g_1 Nmax_6 L_8.npy')

gexact=np.arange(0,2,0.05)

plt.plot(gexact, stark, 'r--')
plt.xlabel(r'$g$')
plt.ylabel(r'$<N_avg>_{GS}$')
red_patch = mpatches.Patch(color='red', label='Exact Diag Wannier stark, h=1, Half filling')

plt.legend(handles=[red_patch] )
plt.show