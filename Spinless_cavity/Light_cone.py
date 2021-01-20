# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 15:18:49 2021

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



def single_el(L):
    ps=['vac']
    for i in range(int(L/2)):
        ps.append('empty')
        ps.append('empty')
    ps[int(L/2)]='full'
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
    
    
    

Nmax, L, g, Omega  =5, 30, 1, 0.5   
dt=0.05
steps=20
J=1
sites = sites(L,Nmax)
ps= product_state(L)
pert=0.2
dt=0.1
psi=MPS.from_product_state(sites, ps)
tmax=3
odd_perm=odd_perm(L)
even_perm=even_perm(L)
delta_t_im=[0.1, 0.01, 0.001, 1.e-4]
chis=[50, 50, 50, 50]
verbose=False
trunc_param=[]
for i in range(4):
    trunc_param.append({'chi_max':chis[i],'svd_min': 0.000000001, 'verbose': verbose})




Peier=npc.outer(npc.expm(-1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']),npc.outer(-J*psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].C.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
Peier_hc=npc.outer(npc.expm(+1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']), npc.outer(-J*psi.sites[1].C.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
cav=npc.outer((Omega/L)*psi.sites[0].N.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
perturbation_ev=npc.outer(psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),pert*psi.sites[1].N.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
perturbation_odd=npc.outer(psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),npc.outer(pert*psi.sites[1].N.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
H_bond=Peier+Peier_hc+cav
H_perturbed_ev=H_bond+perturbation_ev
H_perturbed_odd=H_bond+perturbation_odd


def U_bond(dt, H_bond):
    
    H2 = H_bond.combine_legs([('p0', 'p1', 'p2'), ('p0*', 'p1*', 'p2*')], qconj=[+1, -1])
    H2= (-dt)*H2
    U=npc.expm(H2).split_legs()
    return U




        
def Suz_trot_real_pert(psi, dt, steps):

    U_ev=[U_bond(1j*dt, H_bond)]*int(L/2)
    U_odd= [U_bond((1j*dt)/2, H_bond)]*(int(L/2)-1)
    U_ev[int(L/4)]=U_bond(1j*dt, H_perturbed_ev)
    U_odd[int(L/4)]=U_bond((1j*dt)/2, H_perturbed_odd)
    
    for T in range(steps):
        print("Step number: ", T)

        for i in range(int(L/2)-1): # First Odd sweep
            print("First odd sweep, step: ", i)
            psi.swap_sites(2*i, swap_op=None, trunc_par=trunc_param[1])
            psi.apply_local_op((2*i)+1 , U_odd[i], unitary=True)
            psi.swap_sites(2*i+1, swap_op=None, trunc_par=trunc_param[1])

        for i in range(int(L/2)-1):
            print("ev sweep, step: ", i)

            psi.apply_local_op(L-2-2*i, U_ev[i], unitary=True)
            psi.swap_sites(L-3-2*i, swap_op=None, trunc_par=trunc_param[1])
            psi.swap_sites(L-4-2*i, swap_op=None, trunc_par=trunc_param[1])

        psi.apply_local_op(0, U_ev[0], unitary=True)

        for i in range(int(L/2)-1): # First Odd sweep
            print("Second odd sweep, step: ", i)
            psi.swap_sites(2*i, swap_op=None, trunc_par=trunc_param[1])
            psi.apply_local_op((2*i)+1 , U_odd[i], unitary=True)
            psi.swap_sites(2*i+1, swap_op=None, trunc_par=trunc_param[1])
        for i in range(int(L/2)-1):
            psi.swap_sites(L-3-2*i, swap_op=None, trunc_par=trunc_param[1])
            psi.swap_sites(L-4-2*i, swap_op=None, trunc_par=trunc_param[1])
         
        psi.compress_svd(trunc_param[1])
        print(psi)
        
def Suz_trot_real(psi, dt, steps):
 
    U_ev=U_bond(1j*dt, H_bond)
    U_odd= U_bond((1j*dt)/2, H_bond)
    for T in range(steps):
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


def simple_ev_real(psi, dt, steps):
    U_ev=U_bond(1j*dt, H_bond)
    U_odd= U_bond((1j*dt)/2, H_bond)

    for T in range(steps):
        print("Step number: ", T)

        for i in range(int(L/2)-1): # First Odd sweep
            print("First odd sweep, step: ", i)
            psi.swap_sites(2*i, swap_op=None, trunc_par=trunc_param[1])
            psi.apply_local_op((2*i)+1 , U_odd, unitary=True)
            psi.swap_sites(2*i+1, swap_op=None, trunc_par=trunc_param[1])

        for i in range(int(L/2)-1):
            print("ev sweep, step: ", i)

            psi.apply_local_op(L-2-2*i, U_ev, unitary=True)
            psi.swap_sites(L-3-2*i, swap_op=None, trunc_par=trunc_param[1])
            psi.swap_sites(L-4-2*i, swap_op=None, trunc_par=trunc_param[1])

        psi.apply_local_op(0, U_ev, unitary=True)

        for i in range(int(L/2)-1): # First Odd sweep
            print("Second odd sweep, step: ", i)
            psi.swap_sites(2*i, swap_op=None, trunc_par=trunc_param[1])
            psi.apply_local_op((2*i)+1 , U_odd, unitary=True)
            psi.swap_sites(2*i+1, swap_op=None, trunc_par=trunc_param[1])
        for i in range(int(L/2)-1):
            psi.swap_sites(L-3-2*i, swap_op=None, trunc_par=trunc_param[1])
            psi.swap_sites(L-4-2*i, swap_op=None, trunc_par=trunc_param[1])
         



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
  
    
  

  
    
  
    
  
    
  
    
  
    
  
    
  
with open('psi_ground_L_30g_1J_1.pkl', 'rb') as f:
    psi = pickle.load(f)





n_i_t=[]
n_i=[]
for i in range(L-1):
    n_i.append(psi.expectation_value('N', i+1))
n_i_t.append(n_i)
for i in range(30):
   Suz_trot_real_pert(psi, 0.1, 1)
   n_i=[]
   for j in range(L-1):
       n_i.append(psi.expectation_value('N', j+1))
   n_i_t.append(n_i)
   plt.plot(n_i)
   plt.show()
    
plt.figure()
plt.imshow(n_i_t[::-1],
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(0, L-1, 0, 41))
plt.xlabel('site i')
plt.ylabel('time ')

plt.colorbar().set_label('Occupancy $N$')


"""
delta_t=[0.1]
# Apply a perturbation to PSI
psi.apply_local_op(10, 'N', unitary=True)
psi1=copy.deepcopy(psi)
# Psi1 will be the fixed occupation number to be evolved
psi1.apply_local_op(10, 'N', unitary=True)
for i in range(5):
 Suz_trot_real(psi, delta_t, delta_t)
 Suz_trot_real(psi1, delta_t, delta_t)
 y=[]
 for i in range(L+1):
    psi2=copy.deepcopy(psi)
    y.append(psi1.overlap(psi2.apply_local_op(i, 'N', unitary=True)))
 plt.plot(y)
 plt.ylabel('correlation_function')
 plt.show



"""
