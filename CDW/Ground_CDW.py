# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 08:18:19 2021

@author: giaco
"""

import tenpy
import copy
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/My_library')
from N_cons import Suz_trot_im, H_Peier_bond, product_state, Iterative_g, Iterative_g_from_load
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
from tenpy.linalg.charges import LegCharge, ChargeInfo

def sites(L,Nmax):
 FSite=FermionSite('N', filling=0.5)
 qflat=[[0]]*(Nmax+1)
 ch=ChargeInfo([1], names=None)
 leg = LegCharge.from_qflat(ch, qflat, qconj=1)
 BSite=BosonSite(Nmax=Nmax,conserve=None, filling=0 )
 BSite.change_charge(leg)
 sites=[]
 sites.append(BSite)
 for i in range(L):
     sites.append(FSite)
 return sites

def sites_f(L):
 FSite=FermionSite('N', filling=0.5)

 sites=[]
 
 for i in range(L):
     sites.append(FSite)
 return sites

def product_state_f(L):
    ps=[]
    for i in range(int(L/2)):
        ps.append('empty')
        ps.append('full')
    return ps
J=1
dt=0.05
L=8
g_0=0.1
g=g_0/np.sqrt(L)
Nmax=20
Omega=10
h=0
V=0
mu=0
steps=40
sites_f=sites_f(L)
ps=product_state_f(L)
psi_f=MPS.from_product_state(sites_f, ps)

model_params={'bc_MPS':'finite', 'bc_x':'open', 'explicit_plus_hc':True, 'lattice':'Chain', 'J':J, 'conserve':'N', 'V':V, 'mu':mu, 'L':L}
FC=tenpy.models.fermions_spinless.FermionChain(model_params)

verbose=True
dmrg_params = {
        'mixer': True,  # setting this to True is essential for the 1-site algorithm to work.
        'max_E_err': 1.e-18,
        'trunc_params': {
            'chi_max': 120,
            'svd_min': 1.e-12
        },
        'verbose': verbose,
        'combine': False,
         # specifies single-site
    }
#info = dmrg.run(psi_f, FC, dmrg_params)

with open('Psi_GS_Nmax_20L_8Omega_10J_1h_0V_0g_00.1imTEBD.pkl', 'rb') as f:
    psi = pickle.load(f)

Iterative_g_from_load(psi, 0.1, 1, 0.05, L, Omega, J, h, V, Nmax)


       
     
       
       
       