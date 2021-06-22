# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 16:01:40 2021

@author: giaco
"""

import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/CDW')
import numpy as np
from N_cons import ansatz_wf, Suz_trot_im, H_Peier_bond, ansatz_left, Energy, H_Peier_bond_dn, U_bond, full_sweep, ansatz_first
from tenpy.networks.mps import MPS
import tenpy.linalg.np_conserved as npc
import pickle
import matplotlib.pyplot as plt

Nmax=12
L=24
g_0=0
g=1
Omega  = 10
J=1
h=0
dt=0.001
tmax=5
V=0.0

ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g_0)

V=0.5
trunc_param={'chi_max':80,'svd_min': 1.e-13, 'verbose': False}
psi=ansatz_wf(Nmax, L)


Id=ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
with open('Psi_GS_Nmax_20L_24Omega_10J_1h_0V_0g_00.pkl', 'rb') as f:
    psifermion = pickle.load(f)



for i in range(L):
     psi.set_B(i+1, psifermion.get_B(i))
     psi.set_SL(i+1, psifermion.get_SL(i))
     psi.set_SR(i+1, psifermion.get_SR(i))


   
ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g_0)


occup=[]
for i in range(int(tmax/(10*dt))):
    U=[]
    for z in range(L-1):
      U.append(U_bond(1j*dt, H_Peier_bond_dn(psi, 10*dt*i*g, J/2, Omega,4*V*J, L)))
    print('time_step:', dt*10*i)
    full_sweep(psi, 10, U, Id, trunc_param, L)
    occup.append(psi.expectation_value('N'))
    
np.save('Quench'+ID+'occup.npy', occup)
    


"""
with open(ID+'imTEBD.pkl', 'wb') as f:
       pickle.dump(psi, f)
"""