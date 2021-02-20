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
    for g in [0,0.5,1,2]:
        a=cosa(g,T,Omega)
        b=sum(a)/(100*tmax)
        for i in range(len(a)):
            a[i]=b
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

 
    


Nmax=100
dt=0.01
tmax=5
B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
vac=np.zeros(Nmax+1)
vac[0]=1  
T=-5
g=10
dt=0.01
Omega=10
H=cosm(g*(B+Bd))*T+(Omega*Nb)
U=expm(-1j*dt*H)
X=B+Bd

vac=np.zeros(Nmax+1)
vac[0]=1
v=copy.deepcopy(vac)

def f(x,a,b,c):
    f= a*x+b*np.sqrt(x)+c
    return f 
    
plot_g(0.01, 1, 0.5)
Ts=[]
ggs=[]
for T in -1*np.arange(1,100, 10):
    print(T)
    Ts.append(-T)
    ggs.append(GS(H,Nb, 0.00, 5, 0.01))
plt.plot(Ts, ggs)
a, b = scipy.optimize.curve_fit(f, Ts, ggs)
Ts=np.asarray(Ts)
flist=a[0]*Ts+a[1]*np.sqrt(Ts)+a[2]
plt.plot(Ts, ggs, Ts, flist)








#Plot of cos(A) in quenches were the cavity is populated by states with defined number of photons
"""
CS1=[]
for n in [0,1,5,10,25]:
  gs=np.arange(0,4.05,0.05)
  
  CS0=[]
  v=np.zeros(Nmax+1)
  v[n]=1
  for g in gs:
    CS=cosm(g*(X))
    CS0.append(np.real_if_close(np.tensordot(v.conj().T, np.tensordot(CS, v, 1),1)))
  CS1.append(CS0)
plt.plot(gs,CS1[0],gs,CS1[1],gs,CS1[2],gs,CS1[3],gs,CS1[4])
plt.xlabel(r'$g$')
plt.ylabel(r'$<\cos(A)>_{vac}$')
plt.legend(['n=0','n=1','n=5','n=10','n=25'])
"""

#Plot of cos(A) in quenches were the cavity is populated by states with coherent states
"""
CS1=[]
for alpha in [0,1,5,10,25]:
  gs=np.arange(0,4.01,0.01)
  D=expm(alpha*Bd-np.conj(alpha)*B)
  CS0=[]
  v=np.zeros(Nmax+1)
  v[0]=1
  v=np.tensordot(D, v, 1)
  for g in gs:
    CS=cosm(g*(X))
    CS0.append(np.real_if_close(np.tensordot(v.conj().T, np.tensordot(CS, v, 1),1)))
  CS1.append(CS0)
plt.plot(gs,CS1[0],gs,CS1[1],gs,CS1[2],gs,CS1[3],gs,CS1[4])
plt.xlabel(r'$g$')
plt.ylabel(r'$<\cos(A)>_{vac}$')
plt.legend(['<n>=0','<n>=1','<n>=5','<n>=10','<n>=25'])
"""





        
       

