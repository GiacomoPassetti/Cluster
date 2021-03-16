# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 15:42:21 2021

@author: giaco
"""

import sys 

from ED_functions_fermions_spinless import Vector, states_gen, c_dag_c_L, c_dag_c_i, expectation_value, c_dag_c_i_red, SzSz_L_red, kin_L_red, a_dag_a_red

from scipy.sparse.linalg import eigsh 
import numpy as np
from scipy.sparse.linalg import expm
from time import time

L=10
N=5
N_ph=20
g=2
J=1
U=0.8
Omega=float(sys.argv[1])
PBC=1
OBC=0
BC=PBC

ID='Phase_transitionsU_0-5__g_0-2___Omega_'+str(Omega)
NNN=a_dag_a_red(N_ph, L, N).dot(a_dag_a_red(N_ph, L, N))

Us=np.arange(0,5.2, 0.1)
gg=np.arange(0,2, 0.1)
deviation=np.zeros((len(Us),len(gg)))
n_ph=np.zeros((len(Us),len(gg)))
n_ph_sq=np.zeros((len(Us),len(gg)))

for j in range(len(Us)):
    for k in range(len(gg)):

        
        print('oki', j , k)
        H=-(J/2)*kin_L_red(BC, L, N, N_ph, gg[k], OBC, PBC)+Us[j]*SzSz_L_red(BC, L, N, N_ph, OBC, PBC)+0.02*c_dag_c_i_red (L, N, N_ph, 4)+Omega*a_dag_a_red(N_ph, L, N)
        w, v= eigsh(H, 1, which='SA')
        gs=Vector(v[:, 0])

        nns=[]
        for i in range(L):
             nns.append((gs.expectation_value(c_dag_c_i_red (L, N, N_ph, i))-0.5)**2)
        deviation[j, k]=sum(nns)
        n_ph[j,k]=gs.expectation_value(a_dag_a_red(N_ph, L, N))
        n_ph_sq[j,k]=gs.expectation_value(NNN)

np.save(ID+'large_deviations___0_5.npy', deviation)
np.save(ID+'n_ph.npy', n_ph)
np.save(ID+'n_ph_sq.npy', n_ph_sq)




