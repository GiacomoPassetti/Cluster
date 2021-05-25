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
from scipy.sparse.linalg import eigsh



  
    



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



def cosa_n(n,L, g, Omega):
    v=np.zeros(Nmax+1)
    v[n]=1
    H=-cosm((g/np.sqrt(L))*X)*(L/np.pi)+(Omega*Nb)
    U=expm(-1j*dt*H)
    cosa=cosm(g/np.sqrt(L)*X)
    return expectation(v,cosa)
    

def cos_avgd(n,L, g, Omega):
    v=np.zeros(Nmax+1)
    v[n]=1
    H=-cosm((g/np.sqrt(L))*X)*(L/np.pi)+(Omega*Nb)
    U=expm(-1j*dt*H)
    cosa=cosm(g/np.sqrt(L)*X)
    cosa_t=[cosa_n(n,L, g, Omega)]
    for i in range(int(tmax/dt)):
        v=U.dot(v)
        cosa_t.append(expectation(v, cosa))
    
    val=np.mean(cosa_t)
    
    return val
 
def cos_driven(L, g, Omega, tmax, dt, nu, F):

    H=-cosm((g/np.sqrt(L))*X)*(L/np.pi)+(Omega*Nb)
    w, v= eigsh(H, 1, which='SA')
    v=v[:,0]
    ts=np.arange(0,tmax, dt)
    a_p_a=[-expectation(v, X)]
    for i in range(len(ts)):
      print("t : ", ts[i])
      H=-cosm((g/np.sqrt(L))*X)*(L/np.pi)+(Omega*Nb)+(F*np.sin(nu*ts[i])*X)
      U=expm(-1j*dt*H)
      v=U.dot(v)
      a_p_a.append(-expectation(v, X))
    
    return a_p_a
    


Nmax=100
dt=0.01

B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
vac=np.zeros(Nmax+1)
X=B+Bd
XX=X.dot(X)

L=600

g=0.1
F=0.3
Omega=50
nu=43
tmax=5




apa=cos_driven(L, g, Omega, tmax, dt, nu, F)
ts=np.arange(0,tmax+dt, dt)
plt.figure(dpi=600)
plt.plot(ts, apa)
    








    
    
    









