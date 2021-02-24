# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 08:18:19 2021

@author: giaco
"""

import tenpy
import copy
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/My_library')
from N_cons import Suz_trot_im, H_Peier_bond, product_state
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
L=20
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
psifermion=MPS.from_product_state(sites_f, ps)

model_params={'bc_MPS':'finite', 'bc_x':'open', 'explicit_plus_hc':True, 'lattice':'Chain', 'J':J, 'conserve':'N', 'V':V, 'mu':mu, 'L':L}
FC=tenpy.models.fermions_spinless.FermionChain(model_params)
print(FC.calc_H_bond()[0])
verbose=True
trunc_param={'svd_min': 0.00000000000001, 'verbose': verbose, 'keys':'sorted'}
options={
            'compression_method': 'SVD',
            'trunc_param': trunc_param,
            'keys':'sorted',
            'verbose': verbose 
            }


ID='GS_J_'+str(J)+'V_'+str(V)+'L_'+str(L)


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
info = dmrg.run(psifermion, FC, dmrg_params)

print('Start Imaginary ground search')


plt.plot(psifermion.expectation_value('N'))

print(sites(L, Nmax))
siti = sites(L,Nmax)
ps= product_state(L)




psi=MPS.from_product_state(siti, ps)

for i in range(L):
     psi.set_B(i+1, psifermion.get_B(i))
     psi.set_SL(i+1, psifermion.get_SL(i))
     psi.set_SR(i+1, psifermion.get_SR(i))



max_error_E=[1.e-8, 1.e-7, 1.e-6, 1.e-6, 1.e-6, 1.e-6]
ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g_0)
N_steps=[10, 10, 10, 10, 10, 10]
delta_t_im=[0.1, 1.e-2, 1.e-3, 1.e-4, 1.e-5]
trunc_param={'chi_max':120,'svd_min': 1.e-13, 'verbose': False}

Id=ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
H_bond=[]
for i in range(L-1):
   H_bond.append(H_Peier_bond(psi, g, J, Omega,V, (2*i+1)*h, (2*i+2)*h, L))
   
#Generate the GS from the initial Ansatz
Suz_trot_im(psi, delta_t_im, max_error_E, N_steps, H_bond, trunc_param, L, Id)
plt.plot(psi.expectation_value('N'))
with open(ID+'imTEBD.pkl', 'wb') as f:
       pickle.dump(psi, f)



       
     
       
       
       