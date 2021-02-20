# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 14:28:05 2021

@author: giaco
"""

import tenpy
import copy
import sys
import numpy as np
import numpy.linalg as alg
from tenpy import models
from tenpy.networks.site import FermionSite
from tenpy.networks.site import BosonSite
from tenpy.networks.mps import MPS
from tenpy.tools.params import get_parameter
from tenpy.algorithms.truncation import truncate, svd_theta
import tenpy.linalg.np_conserved as npc
from scipy.linalg import expm
import pickle
import matplotlib.pyplot as plt
import time

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

def mixed_state(L):
    ps=['vac']
    ms = np.array([1/np.sqrt(2), 1/np.sqrt(2)])
    for i in range(int(L/2)):
        ps.append(ms)
        ps.append(ms)
    return ps

def H_Peier(g, J, Omega, U):
    Peier=npc.outer(npc.expm(1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']),npc.outer(-J*psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].C.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
    Peier_hc=npc.outer(npc.expm(-1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']), npc.outer(-J*psi.sites[1].C.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
    cav=npc.outer((Omega/((L-1)))*psi.sites[0].N.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
    rep=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(U*psi.sites[1].N.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].N.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
    H_bond_tebd=Peier+Peier_hc+cav+rep
 #This is the energetic term that will be used in the TEBD algorithm
    
    H_bond=Peier+Peier_hc+rep#This takes into account only the Peier term and is useful to evaluate the energy of the bonds
    return H_bond_tebd, H_bond

def U_bond(dt, H_bond):
    
    H2 = H_bond.combine_legs([('p0', 'p1', 'p2'), ('p0*', 'p1*', 'p2*')], qconj=[+1, -1])
    H2= (-dt)*H2
    U=npc.expm(H2).split_legs()
    return U

def from_full_custom(
                  siti,
                  theta,
                  trunc_par,
                  outer_S,
                  cutoff=1.e-16,
                  form=None,
                  normalize=True,
                  ):

        
        L = len(siti)

        B_list = [None] * L
        S_list = [None] * (L + 1)
        norm = 1. 
        
        labels = ['vL'] + ['p' + str(i) for i in range(L)] + ['vR']
        theta.itranspose(labels)
        # combine legs from left
        for i in range(0, L - 1):
            theta = theta.combine_legs([0, 1])  # combines the legs until `i`
        # now psi has only three legs: ``'(((vL.p0).p1)...p{L-2})', 'p{L-1}', 'vR'``
        for i in range(L - 1, 0, -1):
            # split off B[i]
            theta = theta.combine_legs([labels[i + 1], 'vR'])
            theta, S, B, err, renorm = svd_theta(theta, trunc_par, qtotal_LR=[None, None], inner_labels=['vR', 'vL'])
            
            
            if i > 1:
                theta.iscale_axis(S, 1)
            B_list[i] = B.split_legs(1).replace_label(labels[i + 1], 'p')
            S_list[i] = S
            theta = theta.split_legs(0)
        # psi is now the first `B` in 'A' form
        B_list[0] = theta.replace_label(labels[1], 'p')
        B_form = ['A'] + ['B'] * (L - 1)
        S_list[0], S_list[-1] = outer_S
        res = MPS(siti, B_list, S_list, bc='segment', form=B_form, norm=norm)
        if form is not None:
            res.convert_form(form)
        return res

Nmax=20
L=30
g= 0
Omega  = 10
J=1
h=0
V=0
ps= product_state(L)
sites = sites(L,Nmax)
psi=MPS.from_product_state(sites, ps)
max_error_E=[0.00001, 1.e-5, 1.e-6, 1.e-7, 1.e-8, 1.e-9]
ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)
N_steps=[10, 10, 15, 20, 20, 20]
trunc_param={'chi_max':120,'svd_min': 1.e-13, 'verbose': False}
H_bond_tebd=H_Peier(g, J, Omega, V)[0]
U=U_bond(-0.1, H_bond_tebd)


def apply_local_cav_r(psi, i, op):
            
            cutoff=1.e-13
            "1--  Applico U"

            
            # Prendo il tensore di rango 3 che mi serve e contraggo
            n = 3
            p = psi._get_p_labels(n, False)
            pstar = psi._get_p_labels(n, True)
            th = psi.get_theta(i, n)
            th = npc.tensordot(op, th, axes=[pstar, p])
            
            "2-- Permutazione di indici che realizza lo swap"
            
            th.ireplace_labels(['p0', 'p1','p2'], ['p1', 'p0', 'p2'])
            
            
            
            "3-- Scomposizione in A S B B e ridefinizione di psi"
            split_th = from_full_custom(psi.sites[i:i + n], th, trunc_param,outer_S= (psi.get_SL(i), psi.get_SR(i + n - 1)))
            for j in range(n):
                psi.set_B(i + j, split_th._B[j], split_th.form[j])
            for j in range(n - 1):
                psi.set_SR(i + j, split_th._S[j + 1])
            siteL, siteR = psi.sites[psi._to_valid_index(i)], psi.sites[psi._to_valid_index(i + 1)]
            psi.sites[psi._to_valid_index(i)] = siteR  # swap 'sites' as well
            psi.sites[psi._to_valid_index(i + 1)] = siteL
            


def apply_local_cav_end(psi, i, op):
            cutoff=1.e-13
            "1--  Applico U"

            
            # Prendo il tensore di rango 3 che mi serve e contraggo
            n = 3
            p = psi._get_p_labels(n, False)
            pstar = psi._get_p_labels(n, True)
            th = psi.get_theta(i, n)
            th = npc.tensordot(op, th, axes=[pstar, p])
            

            "3-- Scomposizione in A S B B e ridefinizione di psi"
            split_th = from_full_custom(psi.sites[i:i + n], th, trunc_param,outer_S= (psi.get_SL(i), psi.get_SR(i + n - 1)))
            for j in range(n):
                psi.set_B(i + j, split_th._B[j], split_th.form[j])
            for j in range(n - 1):
                psi.set_SR(i + j, split_th._S[j + 1])


def apply_local_cav_l(psi, i, op):
            i=i-1
            cutoff=1.e-13
            "1--  Genero theta e swappo left"
            n=3
            th = psi.get_theta(i, n)
            th.ireplace_labels(['p0', 'p1','p2'], ['p1', 'p0', 'p2'])
            
            p = psi._get_p_labels(n, False)
            pstar = psi._get_p_labels(n, True)

            
            "2-- applico"

            
            th = npc.tensordot(op, th, axes=[pstar, p])
            
            
            
            
            "3-- Scomposizione in A S B B e ridefinizione di psi"
            split_th = from_full_custom(psi.sites[i:i + n], th, trunc_param,outer_S= (psi.get_SL(i), psi.get_SR(i + n - 1)))
            for j in range(n):
                psi.set_B(i + j, split_th._B[j], split_th.form[j])
            for j in range(n - 1):
                psi.set_SR(i + j, split_th._S[j + 1])
            siteL, siteR = psi.sites[psi._to_valid_index(i)], psi.sites[psi._to_valid_index(i + 1)]
            psi.sites[psi._to_valid_index(i)] = siteR  # swap 'sites' as well
            psi.sites[psi._to_valid_index(i + 1)] = siteL
            




st=time.time()
def full_sweep(step):
  for _ in range(20):
   for i in range(L-1):
     print(i)
     apply_local_cav_r(psi, i, U)   
   apply_local_cav_end(psi, L-2, U)
   for i in range(L-2):
    print(i)
    apply_local_cav_l(psi, L-2-i, U)
   plt.plot(psi.expectation_value('N', list(np.arange(1,L+1))))
   plt.show()
    
  
full_sweep(40)
print(time.time()-st)

