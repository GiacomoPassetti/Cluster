# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 14:31:04 2021

@author: giaco
"""

import tenpy
import copy
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Coherent_state')
from N_cons import sites, product_state, full_sweep, H_Peier_bond, U_bond, single, apply_local_cav_end, full_sweep_second
import numpy as np
import numpy.linalg as alg
import matplotlib.pyplot as plt
from tenpy.networks.site import FermionSite
from tenpy.networks.site import BosonSite
from tenpy.networks.mps import MPS
from tenpy.tools.params import get_parameter
import tenpy.linalg.np_conserved as npc
from scipy.linalg import expm
import pickle
import copy
import time




#Define parameters 
Nmax=2
L=8
g0=0
g= 0/np.sqrt(L)
Omega  = 10
pert=0.1
J=1
h=0
V=0
dt=1/50

tmax=1
N_steps=1
verbose=False
trunc_param={'chi_max':120,'svd_min': 0.00000000000001, 'verbose': verbose}
sites = sites(L,Nmax)
ps1= single(L, 1)
ps2= single(L, 2)

ID='LC_coherent_L'+str(L)+'_g'+str(g)+'_Omega_'+str(Omega)+'dt_'+str(dt)
psi=MPS.from_product_state(sites, ps1)
psi=MPS.from_product_state(sites, ps2)

with open('GS_J_1V_0L_20DMRG.pkl', 'rb') as f:
    psifermion = pickle.load(f)



for i in range(L):
     psi.set_B(i+1, psifermion.get_B(i))
     psi.set_SL(i+1, psifermion.get_SL(i))
     psi.set_SR(i+1, psifermion.get_SR(i))
     



Id=ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
U=[Id]*(L-1)
U[2]=U_bond(-1j*dt/2, H_Peier_bond(psi, g, J, Omega, V, 0, 0, L))
U[3]=U_bond(-1j*dt, H_Peier_bond(psi, g, J, Omega, V, pert, 0, L))

full_sweep(psi, int(tmax/dt), U, Id, trunc_param, L)



print(psi.expectation_value('N')[1:(L+1)])