# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 15:42:21 2021

@author: giaco
"""

import sys 
sys.path.append('C:/Users/giaco/Desktop/Cluster/Shared_Moustafa')
from ED_functions_fermions_spinless import Vector, states_gen, c_dag_c_L, c_dag_c_i, expectation_value, c_dag_c_i_red, SzSz_L_red, kin_L_red, a_dag_a_red
import scipy
from scipy.sparse.linalg import eigsh 
import numpy as np
from scipy.sparse.linalg import expm
from time import time
import matplotlib.pyplot as plt

L=10
N=5
N_ph=15
g=1.75
J=1
U=0
imp=0
Omega=5
PBC=1
OBC=0
BC=PBC
dt=0.1

def H(J, U, g, Omega, imp):
   H=-(J/2)*kin_L_red(BC, L, N, N_ph, g, OBC, PBC)+U*SzSz_L_red(BC, L, N, N_ph, OBC, PBC)+imp*c_dag_c_i_red (L, N, N_ph, 4)+Omega*a_dag_a_red(N_ph, L, N)
   return H
ID='Phase_transitionsU_0-2__g_0-2___Omega_'+str(Omega)
NNN=a_dag_a_red(N_ph, L, N).dot(a_dag_a_red(N_ph, L, N))

H0=H(J, U, g, Omega, imp)
w, v= eigsh(H0, 1, which='SA')
gs=Vector(v[:, 0])

energy=w[0]

Us=np.arange(0.5,1.5,0.01)

e_t=[]
for i in range(len(Us)):
    print('Actual U:', Us[i])
    ID='Phase_transitionsU_0-2__g_0-2___Omega_'+str(Omega)+'U_'+str(Us[i])
    
    H1=H(J, Us[i], g, Omega, imp)
    w, v= eigsh(H1, 1, which='SA')
    H1=H1.toarray()
    U=expm(-1j*dt*H1)
    U=scipy.sparse.csr_matrix(U)
    H1=scipy.sparse.csr_matrix(H1)
    gs.apply(U)
    
    e_t.append(gs.expectation_value(H1)-w[0])
np.save('E_t_slower.npy', e_t)
plt.plot(e_t)
    

    

