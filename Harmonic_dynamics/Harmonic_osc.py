# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 15:28:01 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
from scipy.linalg import expm, sinm, cosm, eigh
from tenpy.networks.site import FermionSite
from tenpy.networks.site import BosonSite
import copy



  
    




    
    





def cosa(g, T, Omega):
    CS=cosm(g*(B+Bd))
    H=cosm(g*(B+Bd))*T+(Omega*Nb)
    U=expm(-1j*dt*H)
    v=copy.deepcopy(vac)
    
    cosa=[np.tensordot(v.conj().T, np.tensordot(CS, v, 1),1)]
    for i in range(int(tmax/dt)):
      v=np.tensordot(U,v,1)
      cosa.append(np.real_if_close(np.tensordot(v.conj().T, np.tensordot(CS, v, 1),1)))
    return cosa
      
def N_av(g, T, Omega):

    H=cosm(g*(B+Bd))*T+(Omega*Nb)
    U=expm(-1j*dt*H)
    v=copy.deepcopy(vac)
    
    cosa=[np.tensordot(v.conj().T, np.tensordot(Nb, v, 1),1)]
    for i in range(int(tmax/dt)):
      v=np.tensordot(U,v,1)
      cosa.append(np.real_if_close(np.tensordot(v.conj().T, np.tensordot(Nb, v, 1),1)))
    return cosa

def X_av(g, T, Omega):

    H=cosm(g*(B+Bd))*T+(Omega*Nb)
    U=expm(-1j*dt*H)
    v=copy.deepcopy(vac)
    X=B+Bd
    cosa=[np.tensordot(v.conj().T, np.tensordot(Nb, v, 1),1)]
    for i in range(int(tmax/dt)):
      v=np.tensordot(U,v,1)
      cosa.append(np.real_if_close(np.tensordot(v.conj().T, np.tensordot(X, v, 1),1)))
    return cosa

def plot_g(gmin, gmax, steps):
    gs=np.arange(gmin, gmax+steps, steps)
    t=np.arange(0,tmax+dt, dt)
    for g in gs:
        a=cosa(g,T,Omega)
        plt.plot(t, a)
    plt.legend(['g=0', 'g=0.5', 'g=1', 'g=1.5', 'g=2'])
    plt.xlabel('t')
    plt.ylabel(r'$<\cos(A(t))>$')
    plt.savefig('Cos(a)_T'+str(T)+'Omega_'+str(Omega)+'.png')
    plt.close('all')

def plot_N(gmin, gmax, steps):
    gs=np.arange(gmin, gmax+steps, steps)
    t=np.arange(0,tmax+dt, dt)
    for g in gs:
        a=N_av(g,T,Omega)
        plt.plot(t, a)
    plt.legend(['g=0', 'g=0.5', 'g=1', 'g=1.5', 'g=2'])
    plt.xlabel('t')
    plt.ylabel(r'$<N(t)>$')
    plt.savefig('Na_T'+str(T)+'Omega_'+str(Omega)+'.png')
    plt.close('all')
    
def plot_X(gmin, gmax, steps):
    gs=np.arange(gmin, gmax+steps, steps)
    t=np.arange(0,tmax+dt, dt)
    for g in gs:
        a=X_av(g,T,Omega)
        plt.plot(t, a)
    plt.legend(['g=0', 'g=0.5', 'g=1', 'g=1.5', 'g=2'])
    plt.xlabel('t')
    plt.ylabel(r'$<X(t)>$')
    plt.savefig('X_T'+str(T)+'Omega_'+str(Omega)+'.png')
    plt.close('all')
    
"""
for T in [-1, -10, -100, -1000]:
 for Omega in [10]:

     plot_X(0, 2.5, 0.5)
"""

Nmax=400
dt=0.01
tmax=3
B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
vac=np.zeros(Nmax+1)
vac[0]=1  
T=-1
g=10
dt=0.01
Omega=10
H=cosm(g*(B+Bd))*T+(Omega*Nb)
U=expm(-1j*dt*H)
X=B+Bd

vac=np.zeros(Nmax+1)
vac[0]=1
v=copy.deepcopy(vac)



gs=np.arange(0,3.01,0.01)
CS0=[]
for g in gs:
    CS=cosm(g*(X))
    CS0.append(np.real_if_close(np.tensordot(v.conj().T, np.tensordot(CS, v, 1),1)))
    
plt.plot(gs,CS0)
plt.xlabel(r'$g$')
plt.ylabel(r'$<Cos(A)>_{vac}$')








        
       

