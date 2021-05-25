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



L=10
hmax= 2
level = 3
J=1
iterations=4
tmax=2
dt=0.01
steps=5
epsilon=0.01
chi=300
k1=1

trunc_param={'chi_max':chi,'svd_min': 1.e-13, 'verbose': False}
options={'N_steps' : steps , 'dt': 0.01, 'trunc_params': trunc_param, 'order': 2 }
ID = "Fidelty_hmax"+str(hmax)+"L_"+str(L)+"iterations_"+str(iterations)+"Level_"+str(level)


       
def ent_state(psi, level, L):
  psi0=copy.deepcopy(psi)
  for j in range(level):
   psi0.apply_local_op(int(L/2), Hadamard)
   psi0.apply_local_op(int(L/2)-1, CNOT)
   for i in range(int(L/2)-1):
    psi0.apply_local_op(int(L/2)-2-i, SWAP)
    psi0.apply_local_op(int(L/2)+i, SWAP)
    
  return psi0

"1) Random Psi"
psi = pol_wf(L, 0)




lattice= Chain(L, psi.sites[0])





ts=np.arange(0, tmax+dt, (dt*steps))
Id, Sx, Sy, Sz, CNOT, SWAP, CZ, Hadamard = Quantum_gates(psi)
ERR_SWAP=epsilon*H_vector(psi, J, 0, 0, L)
ERR_SWAP_sx =epsilon * (npc.outer(Sx.replace_labels(['p', 'p*'], ['p0', 'p0*']),Id.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3]) + npc.outer(Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3]))
sxsx= npc.outer(Sx.replace_labels(['p', 'p*'], ['p0', 'p0*']),Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3])
"Generating the model to perform TEBD"


avg=[]
for j in range(iterations):
   randi=[]
   for z in range(L):
    randi.append(random.random())
   H0=H_bonds_random_randi(psi, 1, L, hmax, randi)
   model=NearestNeighborModel(lattice, H0)
   eng = Engine(psi, model, options)
   psi = pol_wf(L, 1)
   psi1=ent_state(psi, level, L)
   eng1 = Engine(psi1, model, options)
   psi0=copy.deepcopy(psi1)
   fidelty = [abs(psi0.overlap(psi1))]
   for i in range(int(tmax/(steps*dt))):
     eng1.run()
     fidelty.append(abs(psi0.overlap(psi1)))
   avg.append(fidelty)
    
err1 = zip(avg[0], avg[1])

err1 = [x + y for (x, y) in err1] 

for i in range(iterations - 2):
    err1 = zip(err1, avg[2+i])

    err1 = [x + y for (x, y) in err1] 

    
 
err1 = [x/iterations for x in err1]

np.save(ID+"avg_1.npy", err1)
    
    


