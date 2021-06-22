# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 11:35:48 2021

@author: giaco
"""

import sys

sys.path.append('C:/Users/giaco/Desktop/Cluster/Quantum_information')
import numpy as np
from tenpy.networks.mps import MPS
from tenpy.models.model  import NearestNeighborModel
from tenpy.models.lattice import Chain
import tenpy.linalg.np_conserved as npc
import pickle
import random 
from tenpy.algorithms.tebd import Engine
import copy
from tenpy.models.xxz_chain import XXZChain


L = 10
Jx = 1
Delta =  1
steps = 5
tmax = 5
dt = 0.01
model_params = dict(L=L, Jxx=1., Jz=Delta,h = 0,  bc_MPS='finite', conserve='Sz', verbose=0)
chi = 100
trunc_param={'chi_max':chi,'svd_min': 1.e-13, 'verbose': False}
options={'N_steps' : steps , 'dt': 0.01, 'trunc_params': trunc_param, 'order': 2 }
ID = "XXZ_spin_propagation_Delta_"+str(Delta)

M = XXZChain(model_params)
psi = MPS.from_product_state(M.lat.mps_sites(), (["up", "down"] * L)[:L], M.lat.bc_MPS)
eng = Engine(psi, M, options)
a = np.zeros(10)
np.save("provaaa.npy", a)
SZ = []
for i in range(int(tmax/(steps*dt))):
    eng.run()
    SZ.append(psi.expectation_value('Sz'))
np.save(ID+".npy", SZ)



