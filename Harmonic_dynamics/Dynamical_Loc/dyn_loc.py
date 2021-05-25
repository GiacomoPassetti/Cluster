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
 
    


Nmax=400
dt=0.01

B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
vac=np.zeros(Nmax+1)
X=B+Bd
XX=X.dot(X)

L=100
n=20
g=1

Omega=5
tmax=(4*np.pi)/Omega
g_0=g/np.sqrt(L)
nmax=60


cosa=cosm(g/np.sqrt(L)*X)
ns=[1,2,5,10,20,50,100,200]
xs= 2*g_0*np.sqrt(ns)
x=np.arange(0,2*g_0*np.sqrt(ns[-1])+(2*g_0*np.sqrt(ns[-1])/100), 2*g_0*np.sqrt(ns[-1])/100)
y1=[]
for i in ns:
    y1.append(cosa_n(i, L, g, Omega))

Omega=1
y2=[]
for i in ns:
   y2.append(cos_avgd(i, L, g, Omega))

Omega=5
y3=[]
for i in ns:
    y3.append(cos_avgd(i, L, g, Omega))
Omega=10
y4=[]
for i in ns:
    y4.append(cos_avgd(i, L, g, Omega))
y=scipy.special.jv(0,x)
plt.figure(dpi=600)
plt.plot(x, y, "black", xs, y1, "r+", xs, y2, "b.", xs, y3, "g.", xs, y4, "c.")
plt.title("L="+str(L)+"  g="+str(g)+" n=[1,2,5,10,20,50,100,200]")
plt.xlabel("A")
plt.ylabel(r"$W/W_{0}$")

plt.legend([r"$J_{0}(A)$", "<n|Cos(A)|n>", r"$\Omega=1$", r"$\Omega=5$", r"$\Omega=10$"])
    








    
    
    









