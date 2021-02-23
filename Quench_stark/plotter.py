# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:37:14 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt
def datas(L, g_0):
  ID='D:/Data_Spinless_Boson-main/Quench_stark/Quench_ws_Nmax'+str(20)+'L_'+str(L)+'Omega_'+str(10)+'J_'+str(1)+'h_'+str(1)+'V_'+str(0)+'g_0'+str(g_0)
  n=np.load(ID+'n_av.npy')
  nn=np.load(ID+'NN.npy')
  A=np.load(ID+'A.npy')
  eps=np.load(ID+'eps.npy')
  return n, nn, A, eps 


def plot_photons(L, g1, g2, g3):
   n1=datas(L, g1)[1]
   n2=datas(L, g2)[1]
   n3=datas(L, g3)[1]
   t=np.arange(0,10,0.05)
   N1=[]
   for i in range(200):
       N1.append(n1[i][0])
   N2=[]
   for i in range(200):
       N2.append(n2[i][0])
   N3=[]
   for i in range(200):
       N3.append(n3[i][0])

   plt.plot(t,N1, t, N2, t, N3)
   plt.xlabel('t')
   plt.ylabel(r'$<a^{\dagger}a>$')
   plt.legend(['g='+str(g1),'g='+str(g2), 'g='+str(g3) ])
   plt.savefig('Plot_nav_L_'+str(L))
   
n_i_t=datas(30, 2.0)[0]
print(n_i_t[0,30])
n=np.zeros((200,40))
for i in range(200):
    for j in range(40):
        n[i,j]=n_i_t[i,j+1]


        
    
plt.figure()
plt.imshow(n[::-1],
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(0, 40, 0, 10))
plt.xlabel('site i')
plt.ylabel('time g='+str(2))

plt.colorbar().set_label('Occupancy $N$')

   
   