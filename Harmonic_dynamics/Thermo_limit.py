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
import scipy



  
    



def Bessel_0(x):
    val=scipy.special.jv(0,x)
    return val
    
    



def expectation(v, op):
    val = np.real_if_close(np.tensordot(v.conj().T, np.tensordot(op, v, 1),1))
    return val
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



"""
def plot_g(gmin, gmax, steps):
    gs=np.arange(gmin, gmax+steps, steps)
    t=np.arange(0,tmax+dt, dt)
    for g in [0,0.5,1,2]:
        a=cosa(g,T,Omega)
        
        plt.plot(t, a)
    plt.legend(['g=0', 'g=0.5', 'g=1', 'g=2'])
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


def GS(H,Nb, gmin, gmax, steps):
    gs=np.arange(gmin, gmax+steps, steps)
    
    N_ph=[]
    for g in gs:
           H=cosm(g*(B+Bd))*T+(Omega*Nb)
           w, v= eigh(H, eigvals_only=False)
           N_ph.append(np.real_if_close(np.tensordot(v[:,0].conj().T, np.tensordot(Nb, v[:,0], 1),1)))  
    for i in range(len(N_ph)-1):
        
        if N_ph[i+1]<N_ph[i]:
            ggg=gs[i]
            print('break')
            break
    return ggg
"""
 
    


Nmax=200
dt=0.01

B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
vac=np.zeros(Nmax+1)
X=B+Bd

n=20
g=0.5
T=200
Omega=5
tmax=(4*np.pi)/Omega
g_0=g/np.sqrt(T)
nmax=60

H1 = 