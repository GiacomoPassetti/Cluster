# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 13:08:34 2021

@author: giaco
"""

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

def squeezednes_i(g, L):
    H=-cosm((g/np.sqrt(L))*X)*(L/np.pi)+(Omega*Nb)
    w, v= eigsh(H, 1, which='SA')
    v=v[:,0]
    dx=expectation(v, XX)
    dx=dx-(expectation(v, X)**2)

    dy=expectation(v, YY)

    dy=dy-(expectation(v, Y)**2)

    return (dy/dx)
    
   





Nmax=140
L=10000
B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
X= B+Bd
Y= 1j*(Bd-B)
XX=X.dot(X)
YY=Y.dot(Y)
Omega=1
g=np.arange(0,2.025,.025 )
Ls=[10,20,50,100,200,500]
ys1=[]
for L in Ls:

 y1=[]
 for i in range(len(g)):
    y1.append(squeezednes_i(g[i], L))
 ys1.append(y1)
 
Omega=0.1
ys01=[]
for L in Ls:

 y1=[]
 for i in range(len(g)):
    y1.append(squeezednes_i(g[i], L))
 ys01.append(y1)
 
Omega=10
ys10=[]
for L in Ls:

 y1=[]
 for i in range(len(g)):
    y1.append(squeezednes_i(g[i], L))
 ys10.append(y1)
    

    


n_plots=len(ys1)
print(len(ys1))
plt.figure(dpi=600)
# Have a look at the colormaps here and decide which one you'd like:
# http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html



cmap = plt.get_cmap('Greens')
colors = [cmap(i) for i in np.linspace(0.1, 1, n_plots)]

for i, color in enumerate(colors, start=1):
    print(i)
    plt.plot(g, ys1[i-1], color=color)
    
cmap = plt.get_cmap('Blues')
colors = [cmap(i) for i in np.linspace(0.1, 1, n_plots)]

for i, color in enumerate(colors, start=1):
    print(i)
    plt.plot(g, ys01[i-1], color=color)
    
cmap = plt.get_cmap('Oranges')
colors = [cmap(i) for i in np.linspace(0.1, 1, n_plots)]

for i, color in enumerate(colors, start=1):
    print(i)
    plt.plot(g, ys10[i-1], color=color)

plt.xlabel(r"$g$")
plt.ylabel(r"$\frac{\Delta P}{\Delta X}$")
red_patch = mpatches.Patch(color='red', label=r'$\Omega= 10$')
green_patch = mpatches.Patch(color='green', label=r'$\Omega= 1$')
blue_patch = mpatches.Patch(color='blue', label=r'$\Omega= 0.1$')
plt.legend(handles=[blue_patch, green_patch, red_patch])

plt.show()



