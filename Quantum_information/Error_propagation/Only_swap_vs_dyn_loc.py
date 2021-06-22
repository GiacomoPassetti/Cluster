# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 11:10:56 2021

@author: giaco
"""
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quantum_information')
import numpy as np
from functions import random_wf, swap_op, Bond_id, random_swap, H_bonds_random, random_swap_unitary, U_bond, error_op, Quantum_gates, H_bonds_random_randi, H_vector
from tenpy.models.model  import NearestNeighborModel
from tenpy.models.lattice import Chain
import tenpy.linalg.np_conserved as npc
import pickle
import matplotlib.pyplot as plt
import random 
from tenpy.algorithms.tebd import Engine
import copy



L=10

J=0.1

tmax=10
dt=0.01
steps=5
epsilon=0.1
chi=300
k1=1
iterations=10
trunc_param={'chi_max':chi,'svd_min': 1.e-13, 'verbose': False}
options={'N_steps' : steps , 'dt': 0.01, 'trunc_params': trunc_param, 'order': 2, 'tau':0.01 }
ID = "No_err_Swap_J_"+str(J)+"eps_"+str(epsilon)+"Iterations_"+str(iterations)+"t_max"+str(tmax)

       


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

avg1 = []
avg2 = []
avg3 = []
for j in range(iterations):
 psi = random_wf(L)
 randi=[]
 for z in range(L):
    randi.append(random.random())
 H0=H_bonds_random_randi(psi, 0, L, 0, randi)
 H1=H_bonds_random_randi(psi, 0.1, L, 0.02, randi)
 H2=H_bonds_random_randi(psi, 0.1, L, 2, randi)
 H3=H_bonds_random_randi(psi, 0.1, L, 20, randi)
 model0=NearestNeighborModel(lattice, H0)
 model1=NearestNeighborModel(lattice, H1)
 model2=NearestNeighborModel(lattice, H2)
 model3=NearestNeighborModel(lattice, H3)
 print("Iteration n:", j)
 psi0=copy.deepcopy(psi)
 psi1=copy.deepcopy(psi)
 psi2=copy.deepcopy(psi)
 psi3=copy.deepcopy(psi)
 eng0=Engine(psi0, model0, options)
 eng1=Engine(psi1, model1, options)
 eng2=Engine(psi2, model2, options)
 eng3=Engine(psi3, model3, options)
 err1=[psi1.overlap(psi0)]
 err2=[psi2.overlap(psi0)]
 err3=[psi3.overlap(psi0)]

 for i in range(int(tmax/(steps*dt))):
       
       eng0.run()
       eng1.run()
       eng2.run()
       eng3.run()
       
       RI1=random.randint(0, L-2)
       psi0.apply_local_op(RI1, SWAP)
       psi1.apply_local_op(RI1, SWAP + ERR_SWAP)
       psi2.apply_local_op(RI1, SWAP + ERR_SWAP)
       psi3.apply_local_op(RI1, SWAP + ERR_SWAP)

       err1.append(abs(psi0.overlap(psi1)))
       err2.append(abs(psi0.overlap(psi2)))
       err3.append(abs(psi0.overlap(psi3)))
 #plt.plot(ts, err1, ts, err2, ts, err3)
 plt.show()
 avg1.append(err1)
 avg2.append(err2)
 avg3.append(err3)
 
err1 = zip(avg1[0], avg1[1])
err2 = zip(avg2[0], avg2[1])
err3 = zip(avg3[0], avg3[1])
err1 = [x + y for (x, y) in err1] 
err2 = [x + y for (x, y) in err2]
err3 = [x + y for (x, y) in err3]
for i in range(iterations - 2):
    err1 = zip(err1, avg1[2+i])
    err2 = zip(err2, avg2[2+i])
    err3 = zip(err3, avg3[2+i])
    err1 = [x + y for (x, y) in err1] 
    err2 = [x + y for (x, y) in err2]
    err3 = [x + y for (x, y) in err3]
    
 
err1 = [x/iterations for x in err1]
err2 = [x/iterations for x in err2]
err3 = [x/iterations for x in err3]
np.save(ID+"avg_1.npy", err1)
np.save(ID+"avg_2.npy", err2)
np.save(ID+"avg_3.npy", err3)

plt.figure(dpi=800)
ts=np.arange(0, tmax+dt, (dt*steps))
plt.plot( ts, err1, ts, err2, ts, err3)
plt.xlabel("t")
plt.ylabel(r"$<\Psi_{id}|\Psi_{err}>(t)$")
plt.legend([r"$|h_{max}|=0.1J$", r" $|h_{max}|=10J$", r"$|h_{max}|=100J$"])
plt.title(r"Only SWAP  $J=0$  $\epsilon=$"+str(epsilon)+"  Averaged over "+str(iterations) )
#np.save(ID+'SMax.npy', Smax)
#np.save(ID+'Smean.npy', Smean)






