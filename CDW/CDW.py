# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 08:18:15 2021

@author: giaco
"""

import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/CDW')
import numpy as np
from N_cons import ansatz_wf, Suz_trot_im, H_Peier_bond, ansatz_left, Energy, H_Peier_bond_dn

import tenpy.linalg.np_conserved as npc
import pickle
import matplotlib.pyplot as plt


Nmax=12
L=8
g_0=0
g=0
Omega  = 10
J=1/2
h=0

V=2
max_error_E=[1.e-6, 1.e-6, 1.e-3]
ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g_0)
N_steps=[10, 10, 10, 10, 10, 10, 10]
delta_t_im=[1.e-1, 1.e-2, 1.e-3]
trunc_param={'chi_max':120,'svd_min': 1.e-13, 'verbose': False}
psi=ansatz_wf(Nmax, L)
print(psi.sites[1].dN)

Id=ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])

H_bond=[]
for i in range(L-1):
   H_bond.append(H_Peier_bond_dn(psi, g, J/2, Omega,4*V*J, L))


#Generate the GS from the initial Ansatz
Suz_trot_im(psi, delta_t_im, max_error_E, N_steps, H_bond, trunc_param, L, Id)


"""
with open(ID+'imTEBD.pkl', 'wb') as f:
       pickle.dump(psi, f)
"""


