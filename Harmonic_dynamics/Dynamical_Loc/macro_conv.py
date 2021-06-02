# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 11:56:22 2021

@author: giaco
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from tenpy.networks.site import BosonSite
from scipy.linalg import expm, sinm, cosm, eigh
from scipy.sparse.linalg import eigsh
import matplotlib.patches as mpatches



def xi(g, Omega):
    xi=-(1/4)*np.log(1+(g**2)/(np.pi*Omega))
    return xi
def expectation(v, op):
    val = np.real_if_close(np.tensordot(v.conj().T, np.tensordot(op, v, 1),1))
    return val

def squeezednes(g, L):
    H=-cosm((g/np.sqrt(L))*X)*(L/np.pi)+(Omega*Nb)
    w, v= eigsh(H, 1, which='SA')
    v=v[:,0]
    dx=expectation(v, XX)
    dx=dx-(expectation(v, X)**2)

    dy=expectation(v, YY)

    dy=dy-(expectation(v, Y)**2)

    return 1-(dx/dy)

def diff_termo(L, g):
    H1=-cosm((g/np.sqrt(L))*X)*(L/np.pi)+(Omega*Nb)
    w, v= eigsh(H1, 1, which='SA')
    v=v[:,0]
    N_ph_all=expectation(v, XX)
    
    diff=N_ph_all


    return diff

def squeezednes_i(g, L):
    H=-cosm((g/np.sqrt(L))*X)*(L/np.pi)+(Omega*Nb)
    w, v= eigsh(H, 1, which='SA')
    v=v[:,0]
    dx=expectation(v, XX)
    dx=dx-(expectation(v, X)**2)

    dy=expectation(v, YY)

    dy=dy-(expectation(v, Y)**2)

    return np.sqrt(dy/dx)

def double_data(g, L):
    H=-cosm((g/np.sqrt(L))*X)*(L/np.pi)+(Omega*Nb)
    w, v= eigsh(H, 1, which='SA')
    v=v[:,0]
    dx=expectation(v, XX)
    dx=dx-(expectation(v, X)**2)

    dy=expectation(v, YY)

    dy=dy-(expectation(v, Y)**2)

    return np.sqrt(dx), np.sqrt(dy)

Nmax=400
L=1000
B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
X= B+Bd
Y= 1j*(Bd-B)
XX=X.dot(X)
YY=Y.dot(Y)
Omega=1
g=10



Ls=np.arange(10,500, 20)
Y=[]
for i in range(len(Ls)):
    Y.append(1/Ls[i])
g=0.5
Y1=[]
for i in range(len(Ls)):
    Y1.append(diff_termo(Ls[i], g))
g=0.0
Y2=[]
for i in range(len(Ls)):
    Y2.append(diff_termo(Ls[i], g))
g=0.5
Y3=[]
for i in range(len(Ls)):
    Y3.append(diff_termo(Ls[i], g))
g=1
Y4=[]
for i in range(len(Ls)):
    Y4.append(diff_termo(Ls[i], g))
g=2
Y5=[]
for i in range(len(Ls)):
    Y5.append(diff_termo(Ls[i], g))
g=5
Y6=[]
for i in range(len(Ls)):
    Y6.append(diff_termo(Ls[i], g))

plt.figure(dpi=800)
plt.plot(Ls, Y2, Ls, Y3, Ls, Y4, Ls, Y5,  Ls, Y6)
plt.xlabel("L")
plt.ylabel(r"$<a^{\dag}a>_{2nd}-<a^{\dag}a>_{Exact}$")
plt.legend(["g=0.2", "g=0.5", "g=1", "g=2", "g=5"])

