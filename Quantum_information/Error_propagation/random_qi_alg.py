# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 11:10:56 2021

@author: giaco
"""
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quantum_information')
import numpy as np
from functions import random_wf, swap_op, Bond_id, random_swap, H_bonds_random, random_swap_unitary, U_bond, error_op
from tenpy.models.model  import NearestNeighborModel
from tenpy.models.lattice import Chain
import tenpy.linalg.np_conserved as npc
import pickle
import matplotlib.pyplot as plt
import random 
from tenpy.algorithms.tebd import Engine
import copy



L=10

J=1
hmax=30
tmax=3
dt=0.01
steps=5
epsilon=0.5
chi=300
trunc_param={'chi_max':chi,'svd_min': 1.e-13, 'verbose': False}
options={'N_steps' : steps , 'dt': 0.01, 'trunc_params': trunc_param, 'order': 2, 'tau':0.01 }


                   


"1) Random Psi"
psi = random_wf(L)

"2) Operator that introduces the error"
op0 = error_op(psi, dt*steps, 0)
op1 = error_op(psi, dt*steps, epsilon)
op2 = error_op(psi, dt*steps, epsilon/10)
op3 = error_op(psi, dt*steps, epsilon/100)


lattice= Chain(L, psi.sites[0])
Hs=H_bonds_random(psi, J, L, hmax)
model=NearestNeighborModel(lattice, Hs)
eng=Engine(psi, model, options)
ID='Mean_entropy_dt'+str(dt)+'t_max'+str(tmax)+'epsilon_'+str(epsilon)+'h_max'+str(hmax)+'chi_max'+str(chi)+'steps_'+str(steps)
Smean=[np.mean(psi.entanglement_entropy())]
Smax=[max(psi.entanglement_entropy())]
err1=[]
err2=[]
err3=[]
for i in range(int(tmax/(steps*dt))):
       
       eng.run()
       psi1=copy.deepcopy(psi)
       psi2=copy.deepcopy(psi)
       psi3=copy.deepcopy(psi)
       
       RI=random.randint(0, L-2)
       swap_op(psi, RI, op0, trunc_param)
       swap_op(psi1, RI, op1, trunc_param)
       swap_op(psi2, RI, op2, trunc_param)
       swap_op(psi3, RI, op3, trunc_param)
       err1.append(np.abs(psi.overlap(psi1)))
       err2.append(np.abs(psi.overlap(psi2)))
       err3.append(np.abs(psi.overlap(psi3)))
       
       Smean.append(np.mean(psi.entanglement_entropy()))
       Smax.append(max(psi.entanglement_entropy()))

ts=np.arange(0, tmax, (dt*steps))
plt.plot( ts, err2, ts, err3)
np.save(ID+'SMax.npy', Smax)
np.save(ID+'Smean.npy', Smean)

print(err3[5])





