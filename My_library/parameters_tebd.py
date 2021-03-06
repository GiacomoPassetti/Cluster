# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 10:57:56 2021

@author: giaco
"""

import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/My_library')
import apply_op_custom as op
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
from tenpy.linalg.charges import LegCharge, ChargeInfo
from tenpy.algorithms.truncation import truncate, svd_theta, TruncationError
import tenpy.linalg.np_conserved as npc
from scipy.linalg import expm
import pickle
import matplotlib.pyplot as plt
import time




Nmax=10
L=30
g_0=2
g=g_0/np.sqrt(L)
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
delta_t_im=[0.1, 1.e-2, 1.e-3, 1.e-4, 1.e-5]
trunc_param={'chi_max':150,'svd_min': 1.e-13, 'verbose': False}
H_bond=[]
for i in range(L-1):
   H_bond.append(H_Peier_bond(g, J, Omega,V, (2*i+1)*h, (2*i+2)*h, L))

Id=npc.outer(psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])