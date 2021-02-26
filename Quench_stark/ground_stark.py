# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 11:10:56 2021

@author: giaco
"""
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quench_stark')
import numpy as np
from N_cons import ansatz_wf, Suz_trot_im, H_Peier_bond, ansatz_left, Energy, Suz_trot_im_second, Suz_trot_im_try
import tenpy.linalg.np_conserved as npc
import pickle
import matplotlib.pyplot as plt


Nmax=12
L=6
g_0=0
g=0
Omega  = 1
J=1
h=1
V=0
max_error_E=[1.e-8, 1.e-7, 1.e-6, 1.e-6, 1.e-6, 1.e-6, 1.e-6]
ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g_0)
N_steps=[10, 10, 10, 10, 10, 10, 10]
delta_t_im=[0.1, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-6]
trunc_param={'chi_max':120,'svd_min': 1.e-13, 'verbose': False}
psi=ansatz_wf(Nmax, L)
plt.plot(psi.expectation_value('N'))
Id=ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
H_bond=[]
for i in range(L-1):
   H_bond.append(H_Peier_bond(psi, g, J, Omega,V, 0,0, L))

#Generate the GS from the initial Ansatz
Suz_trot_im_try(psi, delta_t_im, max_error_E, N_steps, H_bond, trunc_param, L, Id)
plt.plot(psi.expectation_value('N'))
with open(ID+'imTEBD.pkl', 'wb') as f:
       pickle.dump(psi, f)

