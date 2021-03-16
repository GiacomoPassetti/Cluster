# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 15:42:21 2021

@author: giaco
"""

import sys 
sys.path.append('C:/Users/giaco/Desktop/Cluster/Shared_Moustafa')
from ED_functions_fermions_spinless import Vector, states_gen, c_dag_c_L, c_dag_c_i, expectation_value, c_dag_c_i_red, SzSz_L_red, kin_L_red

from scipy.sparse.linalg import eigsh 
import numpy as np
from scipy.sparse.linalg import expm
from time import time

L=10
N=5
N_ph=10
g=0

PBC=1
OBC=0
BC=PBC
A=-0.5*kin_L_red(BC, L, N, N_ph, g, OBC, PBC)


w, v= eigsh(A, 1, which='SA')
gs=Vector(v[:, 0])


tmax=30
dt=0.05
Us=np.arange(0,4+dt,dt/tmax)
print(Us)
diff=[]
energies=[]

ts=time()
for i in Us:
    ID='Quench_U_exact_L'+str(L)+'J_'+str(0.5)+'U_'+str(0-4)+'g_'+str(g)+'N_ph'+str(N_ph)
    print(i, ts-time())
    H=i*SzSz_L_red(BC, L, N, N_ph, OBC, PBC) -0.5*kin_L_red(BC, L, N, N_ph, g, OBC, PBC)
    e_t, ground = eigsh(H, 1, which='SA')
    np.save(ID+'instant_U_'+str(i)+'groundstate.npy', ground)
    U=expm(-1j*dt*H)
    gs.apply(U)
    diff.append(e_t[0]-gs.expectation_value(H))
    energies.append()
    print(e_t[0]-gs.expectation_value(H))
np.save(ID+'diff.npy', diff)
np.save(ID+'energies.npy', energies)
    

    

