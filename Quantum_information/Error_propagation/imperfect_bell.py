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
level = 2
J=1
iterations=4
tmax=2
dt=0.005
steps=5
eps=0.1
chi=300
k1=1

trunc_param={'chi_max':chi,'svd_min': 1.e-13, 'verbose': False}
options={'N_steps' : steps , 'dt': 0.01, 'trunc_params': trunc_param, 'order': 2 }
ID = "Fidelty_hmax"+str(hmax)+"L_"+str(L)+"iterations_"+str(iterations)+"Level_"+str(level)


       
def ent_state_err(psi, level, L, eps):
  Hadamard_err = copy.deepcopy(Hadamard)
  Hadamard_err[0,0]= -(1-eps)/np.sqrt(2*(1+eps**2))
  Hadamard_err[0,1]= (1-eps)/np.sqrt(2*(1+eps**2))
  Hadamard_err[1,0]= (1+eps)/np.sqrt(2*(1+eps**2))
  Hadamard_err[1,1]= (1+eps)/np.sqrt(2*(1+eps**2))
  psi0=copy.deepcopy(psi)
  for j in range(level):
   psi0.apply_local_op(int(L/2), Hadamard_err)
   psi0.apply_local_op(int(L/2)-1, CNOT)
   for i in range(int(L/2)-1):
    psi0.apply_local_op(int(L/2)-2-i, SWAP)
    psi0.apply_local_op(int(L/2)+i, SWAP)
    
  return psi0

def loc_er(psi, L, eps, i):
    Hadamard_err = copy.deepcopy(Hadamard)
    Hadamard_err[0,0]= -(1-eps)/np.sqrt(2*(1+eps**2))
    Hadamard_err[0,1]= (1-eps)/np.sqrt(2*(1+eps**2))
    Hadamard_err[1,0]= (1+eps)/np.sqrt(2*(1+eps**2))
    Hadamard_err[1,1]= (1+eps)/np.sqrt(2*(1+eps**2))
    psi0=copy.deepcopy(psi)
    print(type(psi0))
    psi0.apply_local_op(i, Hadamard_err)
    psi0.apply_local_op(i-1, CNOT)
    
    return psi0



"1) Random Psi"

psi = pol_wf(L, 1)



lattice= Chain(L, psi.sites[0])





ts=np.arange(0, tmax+dt, (dt*steps))
Id, Sx, Sy, Sz, CNOT, SWAP, CZ, Hadamard = Quantum_gates(psi)


"Generating the model to perform TEBD"
psi1 = pol_wf(L, 1)
psi2 = pol_wf(L, 1)


psi1 = loc_er(psi1, L, 1, 5)








randi=[]
for z in range(L):
   randi.append(random.random())
H0=H_bonds_random_randi(psi, 1, L, 0, randi)
model=NearestNeighborModel(lattice, H0)
eng1 = Engine(psi1, model, options)
eng2 = Engine(psi2, model, options)

psi0 = copy.deepcopy(psi1)
print(abs(psi1.overlap(psi2)))
dat = [abs(psi1.overlap(psi2))]
for i in range(200):
    eng1.run()
    eng2.run()
    
    dat.append(abs(psi1.overlap(psi2)))
    


plt.plot(dat)
"""
plt.figure()
plt.imshow(dat[::-1],
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(0, L, 0, tmax))
plt.xlabel('site i')
plt.ylabel('time ')

plt.colorbar().set_label('Sz')
"""