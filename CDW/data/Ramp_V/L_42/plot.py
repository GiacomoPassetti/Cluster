# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 17:47:26 2021

@author: giaco
"""

import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quench_stark/Data')
import numpy as np
from numpy import exp,arange
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import matplotlib.pyplot as plt
from matplotlib import colors


Nmax=12
L=42
g_0=0.0
g=0.0
Omega  = 20
J=1
h=0
V=4.0
dt=0.005
tmax=10
chi=100



ID='C:/Users/giaco/Desktop/Cluster/CDW/data/Ramp_V/L_42/Quench_ramp_V'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'V_'+str(V)+'g_'+str(g)+'dt_'+str(dt)+'chi_max_'+str(chi)
"""
a=np.load(ID+'occup.npy')
ran=range(1,L+1)
occ=np.zeros((int(1/dt),L)) 

print(occ)

for i in range(int(1/dt)):
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
"""

# the function that I'm going to plot



for g in [0.4, 0.8]:
  ID='C:/Users/giaco/Desktop/Cluster/CDW/data/Ramp_V/L_42/Quench_ramp_V'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'V_'+str(V)+'g_'+str(g)+'dt_'+str(dt)+'chi_max_'+str(chi)
  corrs=np.load(ID+'corrs.npy')
  errors=np.load(ID+'error.npy')
  Z=np.zeros((200,42))
  for i in range(200):
    for j in range(42):
      Z[i,j]=corrs[i,0,j+1]
  divnorm=colors.TwoSlopeNorm(vmin=np.amin(Z), vcenter=0., vmax=np.amax(Z))
  plt.figure()
  plt.imshow(Z[::-1],'RdGy',norm=divnorm,
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(1, 42, 0, 10))

  plt.xlabel('site i   g='+str(g))
  plt.ylabel(r'$t$')


  plt.colorbar().set_label('Occupancy $N$')
  plt.savefig(ID+'correlation.png')
  plt.close()
  plt.plot(errors)
  plt.savefig(ID+'error.png')
  plt.close()
  



