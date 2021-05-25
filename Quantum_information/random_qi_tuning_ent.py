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

J=0.01
hmax=50
tmax=5
dt=0.01
steps=5
epsilon=0.5
chi=300
k1=1
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
Hs_perf=H_bonds_random(psi, 0, L, hmax)
model=NearestNeighborModel(lattice, Hs)
model_perf=NearestNeighborModel(lattice, Hs_perf)
eng=Engine(psi, model, options)
ID='Mean_entropy_dt'+str(dt)+'t_max'+str(tmax)+'epsilon_'+str(epsilon)+'h_max'+str(hmax)+'chi_max'+str(chi)+'steps_'+str(steps)
Smean=[np.mean(psi.entanglement_entropy())]
Smax=[max(psi.entanglement_entropy())]
Stot=[sum(psi.entanglement_entropy())]
err1=[]
err2=[]
err3=[]

Id = psi.sites[0].Id
Sx =2* psi.sites[0].Sx
Sz =2* psi.sites[0].Sz
Sy =2* psi.sites[0].Sy
CNOT = npc.outer(0*psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),0*psi.sites[1].Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3]) 
CNOT[0,0,0,0] = 1
CNOT[1,0,1,0] = 1
CNOT[0,1,1,1] = 1
CNOT[1,1,0,1] = 1
SWAP = npc.outer(0*psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),0*psi.sites[1].Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3]) 
SWAP[1,0,0,1] = 1
SWAP[0,1,1,0] = 1
SWAP[1,1,1,1] = 1
SWAP[0,0,0,0] = 1

CZ = npc.outer(0*psi.sites[0].Id.replace_labels(['p', 'p*'], ['p0', 'p0*']),0*psi.sites[1].Sx.replace_labels(['p', 'p*'], ['p1', 'p1*'])).itranspose([0,2,1,3]) 
CZ[1,0,1,0] = 1
CZ[0,1,0,1] = 1
CZ[0,0,0,0] = 1
CZ[1,1,1,1] = -1

Hadamard =   (1/np.sqrt(2))*(Sx + Sz)  




psi_perf=copy.deepcopy(psi)
qtt=[psi.overlap(psi_perf)]
eng_perf=Engine(psi_perf, model_perf, options)
for i in range(int(tmax/(steps*dt))):
       
       eng.run()
       eng_perf.run()

       
       RI1=random.randint(0, L-2)
       RI2=random.randint(0, L-2)
       RI3=random.randint(0, L-2)
       RI4=random.randint(0, L-2)
       RI5=random.randint(0, L-2)
       psi.apply_local_op(RI2, k1*Sx)
       psi.apply_local_op(RI1, k1*Sy)
       psi.apply_local_op(RI3, k1*CZ)
       psi.apply_local_op(RI4, k1*SWAP)
       psi.apply_local_op(3, Hadamard)

       psi_perf.apply_local_op(RI2, k1*Sx)
       psi_perf.apply_local_op(RI1, k1*Sy)
       psi_perf.apply_local_op(RI3, k1*CZ)
       psi_perf.apply_local_op(RI4, k1*SWAP)
       psi_perf.apply_local_op(3, Hadamard)
       
       Smean.append(np.mean(psi.entanglement_entropy()))
       Smax.append(max(psi.entanglement_entropy()))
       Stot.append(sum(psi.entanglement_entropy()))
       qtt.append(psi.overlap(psi_perf))


ts=np.arange(0, tmax+dt, (dt*steps))
plt.plot( ts, qtt)
#np.save(ID+'SMax.npy', Smax)
#np.save(ID+'Smean.npy', Smean)






