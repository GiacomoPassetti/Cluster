# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 15:42:21 2021

@author: giaco
"""

import sys 
sys.path.append('C:/Users/giaco/Desktop/Cluster/CDW/Quench_Exact')
from ED_functions_fermions_spinless import Vector, states_gen, c_dag_c_L, c_dag_c_i, expectation_value, c_dag_c_i_red, SzSz_L_red, kin_L_red, a_dag_a_red, a_dag_plus_a_red
import scipy
from scipy.sparse.linalg import eigsh 
import numpy as np
from scipy.sparse.linalg import expm
from scipy.linalg import expm as expM
from scipy.linalg import sinm, cosm
from time import time

L=10
N=int(L/2)
N_ph=10
g=0
J=1
U=10
imp=0.1
Omega=10
PBC=1
OBC=0
BC=PBC
dt=0.1

def H(J, U, g, Omega, imp):
   H=-(J)*kin_L_red(BC, L, N, N_ph, g, OBC, PBC)+U*SzSz_L_red(BC, L, N, N_ph, OBC, PBC)+imp*c_dag_c_i_red (L, N, N_ph, 5)+Omega*a_dag_a_red(N_ph, L, N)
   return H
ID='Phase_transitionsU_0-2__g_0-2___Omega_'+str(Omega)
NNN=a_dag_a_red(N_ph, L, N).dot(a_dag_a_red(N_ph, L, N))
#CS=cosm(a_dag_plus_a_red (N_ph, L, N).toarray())
#SN=sinm(a_dag_plus_a_red (N_ph, L, N).toarray())
H0=H(J, U, g, Omega, imp)
w, v= eigsh(H0, 1, which='SA')
gs=Vector(v[:, 0])

energy=w[0]

nns=[]
ns=[]
for i in range(L):
             nns.append((gs.expectation_value(c_dag_c_i_red (L, N, N_ph, i))-0.5)**2)
             ns.append(gs.expectation_value(c_dag_c_i_red (L, N, N_ph, i)))
print(sum(nns))





    

