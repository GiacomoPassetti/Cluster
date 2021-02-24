# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 17:00:40 2021

@author: giaco
"""

import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quench_stark')
import numpy as np
from N_cons import ansatz_wf, Suz_trot_im, H_Peier_bond, full_sweep, U_bond, ansatz_last, ansatz_first, ansatz_left
import tenpy.linalg.np_conserved as npc
import pickle
import time
import matplotlib.pyplot as plt

Nmax=10
L=8
g_0=0
g=g_0/np.sqrt(L)
Omega  = 10
J=1
h=5
V=0
dt=0.005
tmax=10
ID_gs='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(0.1)+'V_'+str(V)+'g_0'+str(g_0)+'imTEBD'
ID='Quench_ws_Nmax'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g_0)

trunc_param={'chi_max':120,'svd_min': 1.e-13, 'verbose': False}
psi=ansatz_left(Nmax, L)
Id=ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
U=[]
for i in range(L-1):
   U.append(U_bond(-1j*dt, H_Peier_bond(psi, g, J, Omega,V, h*(2*i+1), h*(2*i+2), L)))

#Here i will have to load the psi from the ground_data


ts=time.time()
n_av=[]
NN=[]
A=[]
eps=0
errors=[]
n_av.append(psi.expectation_value('N'))
NN.append(psi.expectation_value('NN', [0]))
A.append(psi.expectation_value(psi.sites[0].B+psi.sites[0].Bd, [0]))
errors.append(eps)
for i in range(int(tmax//(10*dt))):
    print(psi.expectation_value('N', ))
    eps += full_sweep(psi, 10, U, Id, trunc_param, L).eps
    plt.plot(psi.expectation_value('N',[1,2,3,4,5,6,7,8]))
    plt.show()

    
print(time.time()-ts)