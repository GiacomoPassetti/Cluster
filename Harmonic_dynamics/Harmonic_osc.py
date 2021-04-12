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

def cosa_n(n,T, g, Omega):
    v=np.zeros(Nmax+1)
    v[n]=1
    H=-cosm((g/np.sqrt(T))*X)*T+(Omega*Nb)
    U=expm(-1j*dt*H)
    cosa=cosm(g/np.sqrt(T)*X)
    return expectation(v,cosa)
    

def cos_avgd(n,T, g, Omega):
    v=np.zeros(Nmax+1)
    v[n]=1
    H=-cosm((g/np.sqrt(T))*X)*T+(Omega*Nb)
    U=expm(-1j*dt*H)
    cosa=cosm(g/np.sqrt(T)*X)
    cosa_t=[cosa_n(n,T, g, Omega)]
    for i in range(int(tmax/dt)):
        v=U.dot(v)
        cosa_t.append(expectation(v, cosa))
    
    val=np.mean(cosa_t)
    
    return val


    

ts=np.arange(0,tmax+dt, dt)

"""
for i in range(int(tmax/dt)):
    v=U.dot(v)
    cosa_t.append(expectation(v, cosa))

plt.plot(ts, cosa_t)   
"""
v=np.zeros(Nmax+1)
v[n]=1
H=-cosm((g/np.sqrt(T))*X)*T+(Omega*Nb)
U=expm(-1j*dt*H)
cosa=cosm(g/np.sqrt(T)*X)
ns=np.arange(0,nmax+1,1)
xs= 2*g_0*np.sqrt(ns)
x=np.arange(0,2*g_0*np.sqrt(nmax)+(2*g_0*np.sqrt(nmax)/100)+(2*g_0*np.sqrt(nmax)/100), 2*g_0*np.sqrt(nmax)/100)
y1=[]
for i in ns:
    y1.append(cosa_n(i, T, g, Omega))

Omega=1
y2=[]
for i in ns:
   y2.append(cos_avgd(i, T, g, Omega))

Omega=5
y3=[]
for i in ns:
    y3.append(cos_avgd(i, T, g, Omega))
Omega=10
y4=[]
for i in ns:
    y4.append(cos_avgd(i, T, g, Omega))
y=scipy.special.jv(0,x)
plt.figure(dpi=600)
plt.plot(x, y, "black", xs, y1, "r+", xs, y2, "b.", xs, y3, "g.", xs, y4, "c.")
plt.title("T="+str(T)+"  g="+str(g))

plt.legend([r"$J_{0}(2 \frac{g}{\sqrt{T} }\sqrt{n})$", "<n|Cos(A)|n>", r"$\Omega=1$", r"$\Omega=5$", r"$\Omega=10$"])
    




"""
Ts=[10, 20, 50]
NN=[]
for i in Ts: 
 T=i
 g_0=g/np.sqrt(T)
 Os=np.arange(0.05,4,0.05)
 n_ph=[]
 for t in Os:
    Omega=t
    H=-cosm(g_0*(B+Bd))*T+(Omega*Nb)
    w, v = eigh(H)
    gs=v[:,0]
    n_ph.append(expectation(gs, X))
 NN.append(n_ph)
    
plt.figure(dpi=600)
plt.plot(Os, NN[0], Os, NN[1], Os, NN[2])
plt.xlabel(r"$\Omega$")
plt.ylabel(r"$<a^{\dagger}a>$")
plt.title(r"$g=$"+str(g)+r"$/ \sqrt{L} $"+"  L="+str(T))
"""

"""
Ts=np.arange(10,5000,20)


gs=[0.2, 0.5, 1, 2]
Omega=5
NN=[]
for i in gs:
 g=i   
 n_ph=[]
 for t in Ts:
    
    g_0=g/np.sqrt(t)
    H=-cosm(g_0*(B+Bd))*t+(Omega*Nb)+sinm(g_0*(B+Bd))
    w, v = eigh(H)
    gs=v[:,0]
  
    n_ph.append(expectation(gs, Nb))
 NN.append(n_ph)
 
plt.figure(dpi=600)
plt.plot(Ts, NN[0], Ts, NN[1], Ts, NN[2], Ts,  NN[3])
plt.xlabel("L")
plt.ylabel(r"$<a^{\dagger}+a>$")
plt.legend(["g=0.2", "g=0.5", "g=1", "g=2"])
plt.title(r"$\Omega=$"+str(Omega)+" J=0.01T ")
"""




    
    
    








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



""" BESSEL"""
"""
xs=np.arange(0,20,0.1)
y=scipy.special.jv(0,xs)
plt.figure(dpi=600)
plt.plot(xs, y, "black")
"""     

