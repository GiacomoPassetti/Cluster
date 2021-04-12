# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 11:10:56 2021

@author: giaco
"""
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quantum_information')
import numpy as np
from functions import random_wf, H_spin, swap_op, Bond_id, random_swap, H_bonds_random, random_swap_unitary
from tenpy.models.model  import NearestNeighborModel
from tenpy.models.lattice import Chain
import tenpy.linalg.np_conserved as npc
import pickle
import matplotlib.pyplot as plt
import random 
from tenpy.algorithms.tebd import Engine



L=40

J=1
hmax=30
tmax=10
dt=0.01
steps=5
epsilon=0.001
chi=300
trunc_param={'chi_max':chi,'svd_min': 1.e-13, 'verbose': False}
options={'N_steps' : steps , 'dt': 0.01, 'trunc_params': trunc_param, 'order': 2, 'tau':0.01 }
psi=random_wf(L)
lattice= Chain(L, psi.sites[0])
Hs=H_bonds_random(psi, J, L, hmax)
model=NearestNeighborModel(lattice, Hs)
eng=Engine(psi, model, options)
ID='Mean_entropy_dt'+str(dt)+'t_max'+str(tmax)+'epsilon_'+str(epsilon)+'h_max'+str(hmax)+'chi_max'+str(chi)+'steps_'+str(steps)
Smean=[np.mean(psi.entanglement_entropy())]
Smax=[max(psi.entanglement_entropy())]
for i in range(int(tmax/(steps*dt))):
       print("aaaaap")
       eng.run()
       random_swap(psi, trunc_param, L, epsilon)
       
       Smean.append(np.mean(psi.entanglement_entropy()))
       Smax.append(max(psi.entanglement_entropy()))
                   
np.save(ID+'SMax.npy', Smax)
np.save(ID+'Smean.npy', Smean)







