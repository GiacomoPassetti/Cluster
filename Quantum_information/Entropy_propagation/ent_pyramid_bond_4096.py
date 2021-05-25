# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 11:10:56 2021

@author: giaco
"""
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quantum_information')
import numpy as np
from tenpy.networks.mps import MPS
from functions import random_wf, swap_op, Bond_id, random_swap, H_bonds_random, random_swap_unitary, U_bond, error_op, Quantum_gates, H_bonds_random_randi, H_vector, pol_wf
from tenpy.models.model  import NearestNeighborModel
from tenpy.models.lattice import Chain
import tenpy.linalg.np_conserved as npc
import pickle
import matplotlib.pyplot as plt
import random 
from tenpy.algorithms.tebd import Engine
import copy



L=50

J=0.01

tmax=10
dt=0.01
steps=5
epsilon=0.01
chi=300
k1=1
iterations=20
trunc_param={'chi_max':chi,'svd_min': 1.e-13, 'verbose': False}
options={'N_steps' : steps , 'dt': 0.01, 'trunc_params': trunc_param, 'order': 2 }
ID = "Random_algorithm_J_Error_sx_"+str(J)+"eps_"+str(epsilon)+"Iterations_"+str(iterations)+"t_max"+str(tmax)

       


"1) Random Psi"
psi = pol_wf(L, 0)




lattice= Chain(L, psi.sites[0])





ts=np.arange(0, tmax+dt, (dt*steps))
Id, Sx, Sy, Sz, CNOT, SWAP, CZ, Hadamard = Quantum_gates(psi)
ERR_SWAP=epsilon*H_vector(psi, J, 0, 0, L)
ERR_SWAP_sx =epsilon * (npc.outer(Sx.replace_labels(['p', 'p*'], ['p0', 'p0*']),Id.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3]) + npc.outer(Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3]))
sxsx= npc.outer(Sx.replace_labels(['p', 'p*'], ['p0', 'p0*']),Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
"Generating the model to perform TEBD"
randi=[]
for z in range(L):
    randi.append(random.random())
H0=H_bonds_random_randi(psi, 0, L, 0, randi)
model=NearestNeighborModel(lattice, H0)
eng = Engine(psi, model, options)

psi = pol_wf(L, 1)

for j in range(3):
  psi.apply_local_op(25, Hadamard)
  psi.apply_local_op(24, CNOT)
  for i in range(22-(2*j)):
    psi.apply_local_op(23-i, SWAP)
    psi.apply_local_op(25+i, SWAP)
    print(psi)






plt.figure(dpi=800)
plt.plot(psi.entanglement_entropy())
plt.xlabel("site i")
plt.ylabel(r"$S_{ab}$")









