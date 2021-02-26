# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 07:45:32 2021

@author: giaco
"""

import tenpy
import copy
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/My_library')
from N_cons import sites, product_state, full_sweep, H_Peier_bond, U_bond
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
Nmax=8
L=20
g0=1
g= 1/np.sqrt(L)
Omega  = 10
pert=0.05
J=2   
h=0
V=0

alpha=0
dt=1/50
tmax=2
N_steps=1
verbose=False
trunc_param={'chi_max':50,'svd_min': 0.00000000000001, 'verbose': verbose}
sites = sites(L,Nmax)
ps= product_state(L)
ID='LC_coherent_L'+str(L)+'_g'+str(g)+'_Omega_'+str(Omega)+'displacement_'+str(alpha)+'dt_'+str(dt)

with open('GS_J_1V_0L_20DMRG.pkl', 'rb') as f:
    psifermion = pickle.load(f)

psi=MPS.from_product_state(sites, ps)

for i in range(L):
     psi.set_B(i+1, psifermion.get_B(i))
     psi.set_SL(i+1, psifermion.get_SL(i))
     psi.set_SR(i+1, psifermion.get_SR(i))

Id=ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
U=[]
for i in range(L-1):
   U.append(U_bond(-1j*dt, H_Peier_bond(psi, g, J, Omega,V, (2*i+1)*h, (2*i+2)*h, L)))
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
    eps += full_sweep(psi, 10, U, Id, trunc_param, L).eps
    plt.plot(psi.expectation_value('N'))
    plt.show()
    n_av.append(psi.expectation_value('N'))
    NN.append(psi.expectation_value('NN', [0]))
    A.append(psi.expectation_value(psi.sites[0].B+psi.sites[0].Bd, [0]))
    errors.append(eps)
    
np.save(ID+'n_av.npy', n_av)
np.save(ID+'NN.npy', NN)
np.save(ID+'A.npy', A)
np.save(ID+'eps.npy', errors)
    
print(time.time()-ts)
