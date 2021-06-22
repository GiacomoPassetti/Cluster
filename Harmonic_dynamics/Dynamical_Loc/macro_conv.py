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
    N_ph_all=expectation(v, Nb)
    
    diff=N_ph_all


    return diff

def diff_termo_2(L, g):
    H1=((g**(2))/(2*np.pi))*XX+(Omega*Nb)
    w, v= eigsh(H1, 1, which='SA')
    v=v[:,0]
    N_ph_all=expectation(v, Nb)
    
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



Ls=np.arange(10,2000, 20)
Y=[]

g=0.2
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
g=3
Y6=[]
for i in range(len(Ls)):
    Y6.append(diff_termo(Ls[i], g))
g=4
Y7=[]
for i in range(len(Ls)):
    Y7.append(diff_termo(Ls[i], g))
g=5
Y8=[]
for i in range(len(Ls)):
    Y8.append(diff_termo(Ls[i], g))


g=0.2
YY2=[]
for i in range(len(Ls)):
    YY2.append(diff_termo_2(Ls[i], g))
g=0.5
YY3=[]
for i in range(len(Ls)):
    YY3.append(diff_termo_2(Ls[i], g))
g=1
YY4=[]
for i in range(len(Ls)):
    YY4.append(diff_termo_2(Ls[i], g))
g=2
YY5=[]
for i in range(len(Ls)):
    YY5.append(diff_termo_2(Ls[i], g))
g=3
YY6=[]
for i in range(len(Ls)):
    YY6.append(diff_termo_2(Ls[i], g))
g=4
YY7=[]
for i in range(len(Ls)):
    YY7.append(diff_termo_2(Ls[i], g))
g=5
YY8=[]
for i in range(len(Ls)):
    YY8.append(diff_termo_2(Ls[i], g))

colors = plt.cm.terrain(np.linspace(0,1,10))
plt.figure(dpi=1200)
plt.plot(Ls, Y2, color = colors[0], label = "g=0.2", linewidth = 2)
plt.plot(Ls, YY2, color = 'black', linestyle = ':', linewidth = 1)
plt.plot(Ls, Y3, color = colors[1], label = "g=0.5", linewidth = 2)
plt.plot(Ls, YY3, color = 'black', linestyle = ':', linewidth = 1)
plt.plot(Ls, Y4, color = colors[2], label = "g=1", linewidth = 2)
plt.plot(Ls, YY4, color = 'black', linestyle = ':', linewidth = 1)
plt.plot(Ls, Y5, color = colors[3], label = "g=2", linewidth = 2)
plt.plot(Ls, YY5, color = 'black', linestyle = ':', linewidth = 1)
plt.plot(Ls, Y6, color = colors[4], label = "g=3", linewidth = 2)
plt.plot(Ls, YY6, color = 'black', linestyle = ':', linewidth = 1)
plt.plot(Ls, Y7, color = colors[5], label = "g=4", linewidth = 2)
plt.plot(Ls, YY7, color = 'black', linestyle = ':', linewidth = 1)
plt.plot(Ls, Y8, color = colors[6], label = "g=5", linewidth = 2)
plt.plot(Ls, YY8, color = 'black', linestyle = ':', linewidth = 1)

plt.xlabel("L")
plt.ylabel(r"$<a^{\dag}a>$")
plt.legend()



