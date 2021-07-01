# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 11:10:56 2021

@author: giaco
"""

import numpy as np
from functions import random_wf, swap_op, Bond_id, random_swap, H_bonds_random, random_swap_unitary, U_bond, error_op, Quantum_gates, H_bonds_random_randi, H_vector
from tenpy.models.model  import NearestNeighborModel
from tenpy.models.lattice import Chain
import tenpy.linalg.np_conserved as npc
import pickle
import sys
import random 
from tenpy.algorithms.tebd import Engine
import copy



L=100

J=1
hmax = float(sys.argv[3])
tmax=5
dt=0.01
steps=5
epsilon=float(sys.argv[2])
chi=200
k1=1
iterations=5
IT = int(sys.argv[1])
trunc_param={'chi_max':chi,'svd_min': 1.e-13, 'verbose': False}
options={'N_steps' : steps , 'dt': 0.01, 'trunc_params': trunc_param, 'order': 2, 'tau':0.01 }
ID = "Random_algorithm_J_Error_sx_"+str(J)+"eps_"+str(epsilon)+"Iterations_"+str(iterations)+"t_max"+str(tmax)+"Random_IT_"+str(IT)+"hMAX_"+str(hmax)

       


"1) Random Psi"
psi = random_wf(L)

"2) Operator that introduces the error"
op0 = error_op(psi, dt*steps, 0)
op1 = error_op(psi, dt*steps, epsilon)
op2 = error_op(psi, dt*steps, epsilon/10)
op3 = error_op(psi, dt*steps, epsilon/100)


lattice= Chain(L, psi.sites[0])





ts=np.arange(0, tmax+dt, (dt*steps))
Id, Sx, Sy, Sz, CNOT, SWAP, CZ, Hadamard = Quantum_gates(psi)
ERR_SWAP=epsilon*H_vector(psi, J, 0, 0, L)
ERR_SWAP_sx =epsilon * (npc.outer(Sx.replace_labels(['p', 'p*'], ['p0', 'p0*']),Id.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3]) + npc.outer(Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3]))
avg1 = []
avg2 = []
avg3 = []
for j in range(iterations):
 psi = random_wf(L)
 randi=[]
 for z in range(L):
    randi.append(random.random())
 H0=H_bonds_random_randi(psi, J, L, 0.02, randi)
 H1=H_bonds_random_randi(psi, J, L, 0.02, randi)

 model0=NearestNeighborModel(lattice, H0)
 model1=NearestNeighborModel(lattice, H1)

 print("Iteration n:", j)
 psi0=copy.deepcopy(psi)
 psi1=copy.deepcopy(psi)
 eng0=Engine(psi0, model0, options)
 eng1=Engine(psi1, model1, options)
 
 err1=[psi1.overlap(psi0)]


 for i in range(int(tmax/(steps*dt))):
       
       eng0.run()
       eng1.run()

       
       RI1=random.randint(0, L-2)
       psi0.apply_local_op(RI1, SWAP)
       psi1.apply_local_op(RI1, SWAP + ERR_SWAP_sx)


       err1.append(abs(psi0.overlap(psi1)))


 avg1.append(err1)

 
err1 = zip(avg1[0], avg1[1])

err1 = [x + y for (x, y) in err1] 

for i in range(iterations - 2):
    err1 = zip(err1, avg1[2+i])

    err1 = [x + y for (x, y) in err1] 

    
 
err1 = [x/iterations for x in err1]

np.save(ID+"avg_1.npy", err1)









