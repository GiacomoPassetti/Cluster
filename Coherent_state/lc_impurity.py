# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 17:12:12 2021

@author: giaco
"""

import tenpy
import copy
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Coherent_state')
from N_cons import sites, product_state, full_sweep, H_Peier_bond, U_bond, single
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
Nmax=1
L=8
g0=0
g= 1/np.sqrt(L)
Omega  = 10
pert=0
J=1
h=0
V=0


dt=1/200
tmax=2
N_steps=1
verbose=False
trunc_param={'chi_max':120,'svd_min': 0.00000000000001, 'verbose': verbose}
sites = sites(L,Nmax)
ps= single(L)
ID='LC_coherent_L'+str(L)+'_g'+str(g)+'_Omega_'+str(Omega)+'dt_'+str(dt)
psi=MPS.from_product_state(sites, ps)

with open('GS_J_1V_0L_20DMRG.pkl', 'rb') as f:
    psifermion = pickle.load(f)



for i in range(L):
     psi.set_B(i+1, psifermion.get_B(i))
     psi.set_SL(i+1, psifermion.get_SL(i))
     psi.set_SR(i+1, psifermion.get_SR(i))

Id=ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
U=[]
for i in range(L-1):
   U.append(U_bond(-1j*dt, H_Peier_bond(psi, g, J, Omega,V, 0, 0, L)))
   
U[int(L/2)-1]=U_bond(-1j*dt, H_Peier_bond(psi, g, J, Omega,V, 0, pert, L))
U[int(L/2)]=U_bond(-1j*dt, H_Peier_bond(psi, g, J, Omega,V, pert, 0, L))

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
    plt.plot(psi.expectation_value('N',[4,5,6]))
    plt.show()
    print((psi.expectation_value('N', [6])-psi.expectation_value('N', [4])))
    n_av.append(psi.expectation_value('N'))
    NN.append(psi.expectation_value('NN', [0]))
    A.append(psi.expectation_value(psi.sites[0].B+psi.sites[0].Bd, [0]))
    errors.append(eps)
    
np.save(ID+'n_av.npy', n_av)
np.save(ID+'NN.npy', NN)
np.save(ID+'A.npy', A)
np.save(ID+'eps.npy', errors)
    
print(time.time()-ts)


plt.figure()
plt.imshow(n_av[::-1],
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(0, L, 0, tmax))
plt.xlabel('site i')
plt.ylabel('time g='+str(g))

plt.colorbar().set_label('Occupancy $N$')

   

    
    

