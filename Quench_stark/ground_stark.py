# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 11:10:56 2021

@author: giaco
"""
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quench_stark')
import numpy as np
from N_cons import ansatz_wf, Suz_trot_im, H_Peier_bond, ansatz_left, Energy
from gs_stark_exact import GS
import tenpy.linalg.np_conserved as npc
import pickle
import matplotlib.pyplot as plt


Nmax=12
L=8
g_0=0
g=1
Omega  = 10
J=1
h=1
V=0
max_error_E=[1.e-8, 1.e-7, 1.e-6, 1.e-6, 1.e-6, 1.e-6, 1.e-6]
ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g_0)
N_steps=[10, 10, 10, 10, 10, 10, 10]
delta_t_im=[0.1, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-6]
trunc_param={'chi_max':120,'svd_min': 1.e-13, 'verbose': False}
psi=ansatz_left(Nmax, L)
plt.plot(psi.expectation_value('N'))
Id=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
H_bond=[H_Peier_bond(psi, g, J, Omega,V, h,h/2, L)]
for i in range(L-3):
   H_bond.append(H_Peier_bond(psi, g, J, Omega,V, (i+2)*h/2,(i+3)*h/2, L))
H_bond.append(H_Peier_bond(psi, g, J, Omega,V, (L-1)*h/2,L*h, L))
gs, energy=GS(g, J, Omega, h)
n1=gs.Nf()

#Generate the GS from the initial Ansatz
Suz_trot_im(psi, delta_t_im, max_error_E, N_steps, H_bond, trunc_param, L, Id)
with open(ID+'imTEBD.pkl', 'wb') as f:
       pickle.dump(psi, f)
ind=[1,2,3,4,5,6,7,8]
n2=psi.expectation_value('N', [1,2,3,4,5,6,7,8])


plt.plot(ind,n1,ind,n2,'r+')
plt.legend(['Exact', 'TEBD'])
