# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 08:18:19 2021

@author: giaco
"""

import tenpy
import copy
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/CDW')
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
from tenpy.algorithms.tebd import Engine

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

dt=0.05
L=24
g_0=0
g=g_0/np.sqrt(L)
Nmax=20
Omega=10
J=1
h=0
V=0
mu=0
steps=40
sit=sites_f(L)
ps=product_state_f(L)
psi_f=MPS.from_product_state(sit, ps)
ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g_0)
model_params={'bc_MPS':'finite', 'bc_x':'open', 'explicit_plus_hc':True, 'lattice':'Chain', 'J':J/2, 'conserve':'N', 'V':0, 'mu':-float(0), 'L':L}
FC=tenpy.models.fermions_spinless.FermionChain(model_params)


    
    

verbose=True
dmrg_params = {
        'mixer': True,  # setting this to True is essential for the 1-site algorithm to work.
        'max_E_err': 1.e-10,
        'trunc_params': {
            'chi_max': 100,
            'svd_min': 1.e-10
        },
        'combine': False,
        'active_sites': 2  # specifies single-site
    }
info = dmrg.run(psi_f, FC, dmrg_params)
plt.plot(psi_f.expectation_value('dN'))
with open(ID+'.pkl', 'wb') as f:
    pickle.dump(psi_f, f)




       
     
       
       
       