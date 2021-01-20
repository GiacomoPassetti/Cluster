# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 17:12:12 2021

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

def displacement_op(alpha):
    D=npc.expm(alpha*psi.sites[0].Bd-np.conj(alpha)*psi.sites[0].B)
    return D

def mixed_state(L):
    ps=['vac']
    ms = np.array([1/np.sqrt(2), 1/np.sqrt(2)])
    for i in range(int(L/2)):
        ps.append(ms)
        ps.append(ms)
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

def H_Peier(g, J, Omega):
    Peier=npc.outer(npc.expm(1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']),npc.outer(-J*psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].C.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
    Peier_hc=npc.outer(npc.expm(-1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']), npc.outer(-J*psi.sites[1].C.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
    cav=npc.outer((Omega/((L-1)))*psi.sites[0].N.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
   
    H_bond_tebd=Peier+Peier_hc+cav  #This is the energetic term that will be used in the TEBD algorithm
    H_bond=Peier+Peier_hc  

    return H_bond_tebd, H_bond
    




def U_bond(dt, H_bond):
    
    H2 = H_bond.combine_legs([('p0', 'p1', 'p2'), ('p0*', 'p1*', 'p2*')], qconj=[+1, -1])
    H2= (-dt)*H2
    U=npc.expm(H2).split_legs()
    return U  



        
def Suz_trot_real(psi, dt, N_steps, H_bond_tebd):
    trunc_err=tenpy.algorithms.truncation.TruncationError(eps=0.0, ov=1.0)
    U_ev=U_bond(1j*dt, H_bond_tebd)
    U_odd= U_bond((1j*dt)/2, H_bond_tebd)

    
    for T in range(N_steps):
        print("Step number: ", T)

        for i in range(int(L/2)-1): # First Odd sweep
            
            trunc_err += psi.swap_sites(2*i, swap_op=None, trunc_par=trunc_param)
            psi.apply_local_op((2*i)+1 , U_odd, unitary=True)
            trunc_err += psi.swap_sites(2*i+1, swap_op=None, trunc_par=trunc_param)

        for i in range(int(L/2)-1):
           

            psi.apply_local_op(L-2-2*i, U_ev, unitary=True)
            trunc_err += psi.swap_sites(L-3-2*i, swap_op=None, trunc_par=trunc_param)
            trunc_err += psi.swap_sites(L-4-2*i, swap_op=None, trunc_par=trunc_param)

        psi.apply_local_op(0, U_ev, unitary=True)

        for i in range(int(L/2)-1): # First Odd sweep
            
            trunc_err += psi.swap_sites(2*i, swap_op=None, trunc_par=trunc_param)
            psi.apply_local_op((2*i)+1 , U_odd, unitary=True)
            trunc_err += psi.swap_sites(2*i+1, swap_op=None, trunc_par=trunc_param)
        for i in range(int(L/2)-1):
            trunc_err += psi.swap_sites(L-3-2*i, swap_op=None, trunc_par=trunc_param)
            trunc_err += psi.swap_sites(L-4-2*i, swap_op=None, trunc_par=trunc_param)
         
        trunc_err += psi.compress_svd(trunc_param)
        print(psi)
    return trunc_err
        

#Define parameters 
Nmax=8
L=10
g= 2
Omega  = 20
J=1   
alpha=1
dt=0.01
tmax=1
N_steps=1
verbose=False
trunc_param={'chi_max':80,'svd_min': 0.00000000000001, 'verbose': verbose}
sites = sites(L,Nmax)
ps= mixed_state(L)
psi=MPS.from_product_state(sites, ps)

#Generite the lattice and the operator

H_bond_tebd=H_Peier(g, J, Omega)[0]
H_bond=H_Peier(g, J, Omega)[1]
D_a=displacement_op(alpha)
X=psi.sites[0].B+psi.sites[0].Bd
#Initialize a WF with the fermions in the ground and the bosonic site in the Coherent space 

psi.apply_local_op(0, D_a)

#Perform a time evolution and save at every instant the average value of X=B+Bd


x_t=[]
t=np.linspace(0, tmax, int(tmax/dt))
for i in range(int(tmax/dt)):
   print("time_step:", i)
   Suz_trot_real(psi, dt, N_steps, H_bond_tebd)
   x_t.append(psi.expectation_value(X, [0]))
   
plt.plot(t, x_t)
plt.xlabel('t')
plt.ylabel('<X(t)>')
plt.show()











"""
n_i_t=[]
n_i=[]
for i in range(L):
    n_i.append(psi.expectation_value('N', i+1))
n_i_t.append(n_i)
for i in range(int(tmax/dt)):
   print("time_step:", i)
   Suz_trot_real_pert(psi, dt, 1, H_bond_tebd, H_perturbed_ev, H_perturbed_odd)
   n_i=[]
   for j in range(L):
       n_i.append(psi.expectation_value('N', j+1))
   n_i_t.append(n_i)
   plt.plot(n_i)
   plt.show()
    
plt.figure()
plt.imshow(n_i_t[::-1],
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(0, L, 0, 1))
plt.xlabel('site i')
plt.ylabel('time g='+str(g))

plt.colorbar().set_label('Occupancy $N$')

np.save('nit_g'+str(g)+'B_4_single_p'+'Omega_'+str(Omega), n_i_t)
"""