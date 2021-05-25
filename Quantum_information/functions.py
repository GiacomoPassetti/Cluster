# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 20:09:40 2021

@author: giaco
"""
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quench_stark')

import tenpy
import copy
import sys
import numpy as np
import numpy.linalg as alg
from tenpy import models
from tenpy.networks.site import SpinSite
import random
from tenpy.networks.mps import MPS
from tenpy.tools.params import get_parameter
from tenpy.linalg.charges import LegCharge, ChargeInfo
from tenpy.algorithms.truncation import truncate, svd_theta, TruncationError
import tenpy.linalg.np_conserved as npc
from scipy.linalg import expm
import pickle
import matplotlib.pyplot as plt
import time

def sites(L):
 FSite=SpinSite(S=0.5, conserve=None)

 sites=[]

 for i in range(L):
     sites.append(FSite)
 return sites

def Quantum_gates(psi):
    Id = psi.sites[0].Id
    Sx =2* psi.sites[0].Sx
    Sz =2* psi.sites[0].Sz
    Sy =2* psi.sites[0].Sy
    CNOT = npc.outer(0*psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),0*psi.sites[1].Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3]) 
    CNOT[0,0,0,0] = 1
    CNOT[1,0,1,0] = 1
    CNOT[0,1,1,1] = 1
    CNOT[1,1,0,1] = 1
    SWAP = npc.outer(0*psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),0*psi.sites[1].Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3]) 
    SWAP[1,0,0,1] = 1
    SWAP[0,1,1,0] = 1
    SWAP[1,1,1,1] = 1
    SWAP[0,0,0,0] = 1

    CZ = npc.outer(0*psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),0*psi.sites[1].Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3]) 
    CZ[1,0,1,0] = 1
    CZ[0,1,0,1] = 1
    CZ[0,0,0,0] = 1
    CZ[1,1,1,1] = -1

    Hadamard =   (1/np.sqrt(2))*(Sx + Sz)  
    return Id, Sx, Sy, Sz, CNOT, SWAP, CZ, Hadamard



def random_product_state(L):
    ps=[]
    for i in range(L):
        ps.append(random.randint(0,1))
        
    return ps

def Bond_id(psi, L):
    Id=npc.outer(psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
    return Id
def H_lutt(psi, J, h1, h2, L):
    
    
    sp=npc.outer(J*psi.sites[0].Sp.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Sm.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
    sm=npc.outer(J*psi.sites[0].Sm.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Sp.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
    ons_l=npc.outer(h1*psi.sites[0].Sz.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
    ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),h2*psi.sites[1].Sz.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
    H_bond=sp+sm+ons_l+ons_r  
    return  H_bond

def H_vector(psi, J, h1, h2, L):
    
    
    sxx=npc.outer(J*psi.sites[0].Sx.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
    syy=npc.outer(J*psi.sites[0].Sy.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Sy.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
    szz=npc.outer(J*psi.sites[0].Sz.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Sz.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
    ons_l=npc.outer(h1*psi.sites[0].Sz.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
    ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),h2*psi.sites[1].Sz.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
    H_bond=sxx+syy+szz+ons_l+ons_r  
    return  H_bond

def H_bonds_random(psi, J, L, hmax):
    randi=[]
    for i in range(L):
        randi.append((random.random()*hmax)-(hmax/2))
    H_B=[None]
    H_B.append(H_vector(psi, J, randi[0], randi[1]/2, L))
    for i in range(L-3):
        H_B.append(H_vector(psi, J, randi[i+1], randi[i+2]/2, L))
    H_B.append(H_vector(psi, J, randi[L-2], randi[L-1]/2, L))
    return H_B
    
def H_bonds_random_randi(psi, J, L, hmax, randi):
    
    H_B=[None]
    H_B.append(H_vector(psi, J, randi[0]*hmax-hmax/2, (randi[1]*hmax/2)-hmax/2, L))
    for i in range(L-3):
        H_B.append(H_vector(psi, J, randi[i+1], randi[i+2]/2, L))
    H_B.append(H_vector(psi, J, randi[L-2], randi[L-1]/2, L))
    return H_B
        
def error_op(psi, dt, epsilon):
    op = npc.outer(epsilon*psi.sites[0].Sx.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])+npc.outer(epsilon*psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
    op_epsilon=U_bond(1j*dt, op)
    return op_epsilon
    

def U_bond(dt, H_bond):
    
    H2 = H_bond.combine_legs([('p0', 'p1'), ('p0*', 'p1*')], qconj=[+1, -1])
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
        return res, err




def random_wf(L):
    ps= random_product_state(L)
    site= sites(L)
    psi=MPS.from_product_state(site, ps)
    return psi
def pol_wf(L, s):
    ps= [s]*L
    site= sites(L)
    psi=MPS.from_product_state(site, ps)
    return psi

def random_swap(psi, trunc_param, L, epsilon, op):
    inter=random.randint(0, L-2)
    op = Bond_id(psi, L)+npc.outer(epsilon*psi.sites[0].Sx.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])+npc.outer(epsilon*psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
    error=swap_op(psi, inter, op, trunc_param)
    return error

def random_swap_unitary(psi, trunc_param, L, epsilon, dt):
    inter=random.randint(0, L-2)
    op=U_bond(dt, npc.outer(epsilon*psi.sites[0].Sx.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])+npc.outer(epsilon*psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),psi.sites[1].Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3]))
    error=swap_op(psi, inter, op, trunc_param)
    return error


def swap_op(psi, i, op, trunc_param):
            
            cutoff=1.e-13
            "1--  Applico U"

            
            # Prendo il tensore di rango 3 che mi serve e contraggo
            n = 2
            p = psi._get_p_labels(n, False)
            pstar = psi._get_p_labels(n, True)
            th = psi.get_theta(i, n)
            th = npc.tensordot(op, th, axes=[pstar, p])
            
            "2-- Permutazione di indici che realizza lo swap"
            
            th.ireplace_labels(['p0', 'p1'], ['p1', 'p0'])
            
            
            
            "3-- Scomposizione in A S B  e ridefinizione di psi"
            split_th, err = from_full_custom(psi.sites[i:i + n], th, trunc_param,outer_S= (psi.get_SL(i), psi.get_SR(i + n - 1)))
            for j in range(n):
                psi.set_B(i + j, split_th._B[j], split_th.form[j])
            for j in range(n - 1):
                psi.set_SR(i + j, split_th._S[j + 1])
            siteL, siteR = psi.sites[psi._to_valid_index(i)], psi.sites[psi._to_valid_index(i + 1)]
            psi.sites[psi._to_valid_index(i)] = siteR  # swap 'sites' as well
            psi.sites[psi._to_valid_index(i + 1)] = siteL
            
            return err
            












