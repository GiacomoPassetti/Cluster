# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 13:48:53 2021

@author: giaco
"""

import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quench_stark/Data')
import numpy as np

import matplotlib.pyplot as plt

Nmax=15
L=24
g_0=1.0
g=0
Omega  = 10
J=1
h=1
V=0
dt=0.001
tmax=10


"""
c=[]
for g_0 in [0.2,0.3,0.4,0.5]:
    b=[]
    ID='C:/Users/giaco/Desktop/Cluster/Quench_stark/Data_quench/Quench_stark_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g_0)+'dt_'+str(dt)+'tmax'+str(tmax)
    a=np.load(ID+'n_av.npy')
    for i in range(len(a)):
       b.append(a[i][0])
    c.append(b)
    
    

    
    
t=np.arange(0,10,0.01)
gg=[0]+list(np.arange(0.2,2,0.1))
plt.plot(t, c[0],t, c[1],t, c[2],t,c[3])
plt.xlabel(r'$t$')
plt.ylabel(r'$<a^{\dagger}a>$')
plt.legend(['g=0.2','g=0.3','g=0.4','g=0.5'])
"""


"""
    b=[]
    for i in range(len(a)):
       b.append(a[i][0])
plt.plot(ind, b[9])
plt.xlabel(r'$g$')
plt.ylabel(r'$<c^{\dagger}c>$')
plt.legend(['L=24'])
"""
"""
c=[]
for g_0 in [0.2,1.0,1.8]:
    b=[]
    ID='C:/Users/giaco/Desktop/Cluster/Quench_stark/Data_quench/Quench_stark_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g_0)+'dt_'+str(dt)+'tmax'+str(tmax)
    b.append(np.load(ID+'eps.npy').reshape((1000)))
    
    c.append(np.load(ID+'eps.npy'))

t=np.arange(0,10,0.01)
gg=[0]+list(np.arange(0.2,2,0.1))
plt.plot(t, c[0],t, c[1],t, c[2])
plt.xlabel(r'$t$')
plt.ylabel(r'$Error SVD$')
plt.legend(['g=0.2','g=1.0','g=1.8'])
"""


ID='C:/Users/giaco/Desktop/Cluster/Quench_stark/Data_quench/Quench_stark_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g_0)+'dt_'+str(dt)+'tmax'+str(tmax)
a=np.load(ID+'n_av.npy')
ran=range(1,L+1)
occ=np.zeros((1000,L))

for i in range(1000):
    for j in range(L):
        occ[i,j]=a[i][j+1]

plt.figure()
plt.imshow(occ[::-1],'plasma',
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(1, 24, 0, 10))

plt.xlabel('site i   g='+str(g_0))
plt.ylabel(r'$t$')


plt.colorbar().set_label('Occupancy $N$')
