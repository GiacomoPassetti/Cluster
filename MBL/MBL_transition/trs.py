# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 08:48:46 2021

@author: giaco
"""

import sys
sys.path.append('C:/Users/giaco/Desktop/MBL')
import numpy as np
from tenpy.networks.mps import MPS
from functions import random_wf, swap_op, Bond_id, random_swap, H_bonds_random, random_swap_unitary, U_bond, error_op, H_bonds_random_randi, H_vector, pol_wf
from tenpy.models.model  import NearestNeighborModel
from tenpy.models.lattice import Chain
import tenpy.linalg.np_conserved as npc
import pickle
import matplotlib.pyplot as plt
import random 
from tenpy.algorithms.tebd import Engine
import copy


L=30
hmax= 50



tmax=10
dt=0.005
steps=5
epsilon=0.01
chi=150
k1=1
RI = 1

trunc_param={'chi_max':chi,'svd_min': 1.e-13, 'verbose': False}
options={'N_steps' : steps , 'dt': dt, 'trunc_params': trunc_param, 'order': 2 }
ID = "MBLtrs_hmax"+str(hmax)+"L_"+str(L)+"RI_"+str(RI)

psi = random_wf(L)
lattice= Chain(L, psi.sites[0])
ts=np.arange(0, tmax+dt, (dt*steps))


ent=[]


H0=H_bonds_random(psi, 1, L, hmax)
model=NearestNeighborModel(lattice, H0)
eng = Engine(psi, model, options)


ent.append(psi.entanglement_entropy())  
for i in range(int(tmax/(steps*dt))):
     eng.run()
     ent.append(psi.entanglement_entropy())
     
     if max(psi.chi)>120:
         ID = ID + 'ABORTED_AT_t_' + str(i*dt*steps)
         break
     


np.save(ID+'.npy', ent)