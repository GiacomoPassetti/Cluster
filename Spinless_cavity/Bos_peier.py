# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 17:37:26 2021

@author: giaco
"""

import tenpy
import copy
import sys
import numpy as np
import numpy.linalg as alg
import matplotlib.pyplot as plt
from tenpy import models
from tenpy.networks.site import SpinSite
from tenpy.networks.site import FermionSite
from tenpy.networks.site import BosonSite
from tenpy.models.model import CouplingModel
from tenpy.models.model import CouplingMPOModel
from tenpy.models.spins import SpinModel
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS
from tenpy.models.lattice import Lattice
from tenpy.tools.params import get_parameter
import tenpy.linalg.np_conserved as npc
from tenpy.networks.mpo import MPO, MPOEnvironment
import tenpy.linalg.charges as charges
from tenpy.models.lattice import Chain
from scipy.linalg import expm
from tenpy.models.fermions_spinless import FermionModel
from tenpy.algorithms.tebd import Engine
import pickle
import copy

def sites(L,Nmax):
 FSite=FermionSite(None, filling=0.5)
 BSite=BosonSite(Nmax=Nmax,conserve=None, filling=0 )
 sites=[]
 sites.append(BSite)
 for i in range(L):
     sites.append(FSite)
 return sites

def product_state(L):
    ps=['vac']
    for i in range(int(L/2)):
        ps.append('empty')
        ps.append('full')
    return ps

def psi(sites,ps):
    psi=MPS.from_product_state(sites, ps)
    return psi

def odd_perm(L):
    perm=[]
    b=np.arange(L+1)
    b[0]=1
    b[1]=0
    perm.append(b)
    for i in range(int(L/2)-2):
      b=np.arange(L+1)
      b[2*i+1]=2*i+3
      b[2*i+3]=2*i+1
      perm.append(b)
    return perm
def even_perm(L):
    perm=[]
    b=np.arange(L+1)

    perm.append(b)
    for i in range(int(L/2)-1):
      b=np.arange(L+1)
      b[2*i]=2*i+2
      b[2*i+2]=2*i
      perm.append(b)
    return perm
def last_perm(L):
    perm_od=np.arange(L+1)
    perm_ev=np.arange(L+1)
    perm_od[0]=L-3
    perm_od[L-3]=0
    perm_ev[0]=L-2
    perm_ev[L-2]=0
    return perm_od, perm_ev
    
    
    

Nmax, L, g, Omega  =5, 20, 0, 0.5   
dt=0.05
steps=20
J=2
sites = sites(L,Nmax)
ps= product_state(L)
dt=0.1
psi=psi(sites,ps)
tmax=3
odd_perm=odd_perm(L)
even_perm=even_perm(L)
delta_t_im=[0.1, 0.01, 0.001, 1.e-4]
chis=[40, 40, 50, 50]
verbose=False
trunc_param=[]
for i in range(4):
    trunc_param.append({'chi_max':chis[i],'svd_min': 0.000000001, 'verbose': verbose})




Peier=npc.outer(npc.expm(-1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']),npc.outer(-J*psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].C.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
Peier_hc=npc.outer(npc.expm(+1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']), npc.outer(-J*psi.sites[1].C.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
cav=npc.outer(Omega*psi.sites[0].N.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
H_bond=Peier+Peier_hc+cav

def U_bond(dt, H_bond):
    
    H2 = H_bond.combine_legs([('p0', 'p1', 'p2'), ('p0*', 'p1*', 'p2*')], qconj=[+1, -1])
    H2= (-dt)*H2
    U=npc.expm(H2).split_legs()
    return U

def Suz_trot_im(delta_t, max_error_E, N_steps):
 DeltaE=2*max_error_E
 E_old=Energy(psi)[1]
 for dt in range(len(delta_t)):
    print("delta_tau =", delta_t[dt])
    U_ev=U_bond(delta_t[dt], H_bond)
    U_odd= U_bond(delta_t[dt]/2, H_bond)
    DeltaE= 2*max_error_E
    step=0
    while (DeltaE > max_error_E):
      for T in range(N_steps): 

        print("Step:", T)
        for i in range(int(L/2)-1):
            psi.permute_sites(odd_perm[i], swap_op=None, trunc_par=trunc_param[dt])
            psi.apply_local_op((2*i)+1 , U_odd, unitary=False)
        psi.permute_sites(last_perm(L)[0], swap_op=None, trunc_par=trunc_param[dt])
        for i in range(int(L/2)):
            psi.permute_sites(even_perm[i], swap_op=None, trunc_par=trunc_param[dt])
            psi.apply_local_op(2*i, U_ev, unitary=False)
        psi.permute_sites(last_perm(L)[1], swap_op=None, trunc_par=trunc_param[dt])
        for i in range(int(L/2)-1):
            psi.permute_sites(odd_perm[i], swap_op=None, trunc_par=trunc_param[dt])
            psi.apply_local_op((2*i)+1 , U_odd, unitary=False)
        psi.permute_sites(last_perm(L)[0], swap_op=None, trunc_par=trunc_param[dt])
        psi.compress_svd(trunc_param[dt])
      
      step += N_steps
      E=Energy(psi)[1]
      DeltaE=np.abs(E_old-E)
      E_old=E
      
      print("After", step, "steps, bond_E = ", E, "and DeltaE = ", DeltaE )
        
def Suz_trot_real(delta_t, tmax):
 for dt in delta_t:
    U_ev=U_bond(1j*dt, H_bond)
    U_odd= U_bond((1j*dt)/2, H_bond)
    for T in range(int(tmax/dt)):
        """First odd string"""

        print("Step:", int(tmax/dt))
        for i in range(int(L/2)-1):
            psi.permute_sites(odd_perm[i], swap_op=None, trunc_par=trunc_param[1])
            psi.apply_local_op((2*i)+1 , U_odd, unitary=True)
        psi.permute_sites(last_perm(L)[0], swap_op=None, trunc_par=trunc_param[1])
        for i in range(int(L/2)):
            psi.permute_sites(even_perm[i], swap_op=None, trunc_par=trunc_param[1])
            psi.apply_local_op(2*i, U_ev, unitary=True)
        psi.permute_sites(last_perm(L)[1], swap_op=None, trunc_par=trunc_param[1])
        for i in range(int(L/2)-1):
            psi.permute_sites(odd_perm[i], swap_op=None, trunc_par=trunc_param[1])
            psi.apply_local_op((2*i)+1 , U_odd, unitary=True)
        psi.permute_sites(last_perm(L)[0], swap_op=None, trunc_par=trunc_param[1])
        psi.compress_svd(trunc_param[1])
        print(psi)
        
        
def Energy(psi):
        E=[]
        E.append(psi.expectation_value(Omega*sites[0].N, [0]))
        for i in range(int(L/2)-1):
            psi.permute_sites(odd_perm[i], swap_op=None, trunc_par=trunc_param[1])
            E.append(psi.expectation_value(H_bond , [2*i+1]))
        psi.permute_sites(last_perm(L)[0], swap_op=None, trunc_par=trunc_param[1])
        for i in range(int(L/2)):
            psi.permute_sites(even_perm[i], swap_op=None, trunc_par=trunc_param[1])
            E.append(psi.expectation_value(H_bond, [2*i]))
        psi.permute_sites(last_perm(L)[1], swap_op=None, trunc_par=trunc_param[1])
        E_tot=np.sum(E)
        
        E_bond=np.average(E[1:L])
        return E_tot, E_bond
        

    

data = {"psi": psi,  # e.g. an MPS
        
        "parameters": {"L": 30, "g": 0.5, "J":1, "Omega":0.5}}

Suz_trot_im(delta_t_im, 0.00001, 5)

perturb=sites[1].N-(1/2)*sites[1].Id
psi.apply_local_op(15, perturb, unitary=True)


S = [psi.entanglement_entropy()]
for n in range(int(tmax/dt)):
        Suz_trot_real([dt], dt)
        S.append(psi.entanglement_entropy())

plt.figure()
plt.imshow(S[::-1],
               vmin=0.,
               aspect='auto',
               interpolation='nearest',
               extent=(0, L , 0, 3))
plt.xlabel('site $i$')
plt.ylabel('time $t/J$')
plt.ylim(0., tmax)
plt.colorbar().set_label('entropy $S$')
    

