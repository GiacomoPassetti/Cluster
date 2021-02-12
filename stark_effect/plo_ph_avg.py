# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 10:20:47 2021

@author: giaco
"""
import numpy as np
import matplotlib.pyplot as plt 



gsn=list(np.arange(0, 2.2, 0.2))
gs=list(np.arange(0,2,0.2))
n=np.load('C:/Users/giaco/Desktop/Cluster/Onsite/data_gs_ons/Average_boson_half_filling_starkOmega_1J_1 g_1 Nmax_6 L_8.npy')
nt=np.load('C:/Users/giaco/Desktop/Cluster/Onsite/data_gs_ons/N_avg_bos__TEBDPsi_GS_Nmax_8L_8g_1.8Omega_1J_1h_1V_0.npy')
dn=np.zeros(10)
for i in range(10):
    dn[i]=n[i]-(nt[i]**2)

plt.plot(gsn, n, 'r--', gs, nt, 'b')
plt.ylabel(r'$<N>$')
plt.xlabel('g')
plt.legend([r'$ED$', r'$TEBD$'])
plt.show

