# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 11:03:06 2021

@author: giaco
"""

import matplotlib.pyplot as plt
import numpy as np

Omega=10
L=10
a=np.load('C:/Users/giaco/Desktop/Cluster/CDW/Quench_Exact/Phase_transition_data/Phase_transitionsU_0-2__g_0-2___Omega_'+str(Omega)+'.0deviations.npy')
b=np.load('C:/Users/giaco/Desktop/Cluster/CDW/Quench_Exact/Phase_transition_data/Phase_transitionsU_0-2__g_0-2___Omega_10.0n_ph.npy')
print(b.shape, a.shape)


"""
for Omega in [5, 10, 20]:
  ID='C:/Users/giaco/Desktop/Cluster/CDW/Quench_Exact/Phase_transition_data/Phase_transitionsU_0-2__g_0-2___Omega_'+str(Omega)+'.0deviations.npy'
  a=np.load('C:/Users/giaco/Desktop/Cluster/CDW/Quench_Exact/Phase_transition_data/Phase_transitionsU_0-2__g_0-2___Omega_'+str(Omega)+'.0deviations.npy')
  print(a.shape)
  plt.figure()
  plt.imshow(a[::-1],'plasma',
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(0, 2, 0, 2))

  plt.xlabel(r'$U$')
  plt.ylabel(r'$g$')


  plt.colorbar().set_label(r'$\sum_{i}(n_{i}-\frac{1}{2})^{2}$')
  plt.savefig(ID+'plot_'+str(Omega)+'.png')

  plt.close()
"""
  
"""
for Omega in [10]:
  ID='C:/Users/giaco/Desktop/Cluster/CDW/Quench_Exact/Phase_transition_data/Phase_transitionsU_0-2__g_0-2___Omega_'+str(Omega)+'.0deviations.npy'
  a=np.load('C:/Users/giaco/Desktop/Cluster/CDW/Quench_Exact/Phase_transition_data/Phase_transitionsU_0-5__g_0-2___Omega_10.0large_deviations___0_5.npy')
  print(a.shape)
  plt.figure()
  plt.imshow(a[::-1],'plasma',
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(0, 5, 0, 2))

  plt.xlabel(r'$U$')
  plt.ylabel(r'$g$')


  plt.colorbar().set_label(r'$\sum_{i}(n_{i}-\frac{1}{2})^{2}$')
  plt.savefig(ID+'plot_'+str(Omega)+'.png')
  plt.close()
"""

"""
z=np.zeros((41,41))
Us=np.arange(0,4.1,0.1)
gg=np.arange(0,2.05,0.05)
for j in range(41):
    for k in range(41):

        z[j,k]=Us[j]*gg[k]
        
plt.figure()
plt.imshow(z[::-1],'plasma',
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(0, 2, 0, 2))

plt.xlabel(r'$g$')
plt.ylabel(r'$U$')


plt.colorbar().set_label(r'$\sum_{i}(n_{i}-\frac{1}{2})^{2}$')
"""

"""
ID='C:/Users/giaco/Desktop/Cluster/CDW/Quench_Exact/Phase_transition_data/Phase_transitionsU_0-2__g_0-2___Omega_'+str(Omega)+'.0deviations.npy'
a=np.load('C:/Users/giaco/Desktop/Cluster/CDW/Quench_Exact/Phase_transition_data/Phase_transitionsU_0-5__g_0-2___Omega_10.0L_12cd large_deviations___0_5.npy')
print(a[0,0])
plt.figure()
plt.imshow(a[::-1],'plasma',
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(0, 2, 0, 5))

plt.xlabel(r'$g$')
plt.ylabel(r'$U$')


plt.colorbar().set_label(r'$\sum_{i}(n_{i}-\frac{1}{2})^{2}$')
plt.savefig(ID+'final_plot_Omega_'+str(Omega)+'L_'+str(L)+'.png')
"""

