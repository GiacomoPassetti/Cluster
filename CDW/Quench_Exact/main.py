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
g=0.5
J=2
U=0
imp=0.0
Omega=10
PBC=1
OBC=0
BC=PBC
dt=0.1

def H(J, U, g, Omega, imp):
   H=-(J/2)*kin_L_red(BC, L, N, N_ph, g, OBC, PBC)+U*SzSz_L_red(BC, L, N, N_ph, OBC, PBC)+imp*c_dag_c_i_red (L, N, N_ph, 5)+Omega*a_dag_a_red(N_ph, L, N)
   return H
ID='Phase_transitionsU_0-2__g_0-2___Omega_'+str(Omega)
NNN=a_dag_a_red(N_ph, L, N).dot(a_dag_a_red(N_ph, L, N))
CS=cosm(a_dag_plus_a_red (N_ph, L, N).toarray())
SN=sinm(a_dag_plus_a_red (N_ph, L, N).toarray())
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
print(ns)
print('Cos(A)=', gs.expectation_value(CS))

print('Sin(A)=', gs.expectation_value(SN))



"""
print('done', energy)
H1=H(J, U, g, Omega, imp)
H1=H1.toarray()
U=expm(-1j*dt*H1)
U=scipy.sparse.csr_matrix(U)
H1=scipy.sparse.csr_matrix(H1)
print('done')
Us=np.arange(0,6,0.00001)

e_t=[]
for i in range(len(Us)):
    ts=time()
    print('Actual U:', Us[i], 'time:', time()-ts)
    H1=H(J, Us[i], g, Omega, imp)
    w, v= eigsh(H1, 1, which='SA')
    H1=H1.toarray()
    U=expm(-1j*dt*H1)
    U=scipy.sparse.csr_matrix(U)
    H1=scipy.sparse.csr_matrix(H1)
    gs.apply(U)
    
    print(gs.expectation_value(H1)-w[0])
np.save('E_t.npy', e_t)
"""  


"""
Us=np.arange(0,2.05, 0.05)
gg=np.arange(0,2.05, 0.05)
deviation=np.zeros((41,41))
nph=np.zeros((41,41))
nph_sq=np.zeros((41,41))
for U in Us:
    for g in gg:
        j=int(U/0.05)
        k=int( g//0.05)
        print(j,k)
        
        H=-(J/2)*kin_L_red(BC, L, N, N_ph, g, OBC, PBC)+U*SzSz_L_red(BC, L, N, N_ph, OBC, PBC)+0.05*c_dag_c_i_red (L, N, N_ph, 4)+Omega*a_dag_a_red(N_ph, L, N)
        w, v= eigsh(H, 1, which='SA')
        gs=Vector(v[:, 0])
        Us=np.arange(0,2.05, 0.05)
        gg=np.arange(0,2.05, 0.05)
        nns=[]
        for i in range(L):
             nns.append((gs.expectation_value(c_dag_c_i_red (L, N, N_ph, i))-0.5)**2)
        deviation[j, k]=sum(nns)
        nph[j,k]=gs.expectation_value(a_dag_a_red(N_ph, L, N))
        nph_sq[j,k]=gs.expectation_value(NNN)
np.save(ID+'deviations.npy', deviation)
np.save(ID+'n_ph.npy', nph)
"""



#  Quench
"""
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
"""

    

