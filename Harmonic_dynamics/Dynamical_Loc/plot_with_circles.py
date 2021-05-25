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
    
   





Nmax=180
L=10000
B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
X= B+Bd
Y= 1j*(Bd-B)
XX=X.dot(X)
YY=Y.dot(Y)
Omega=1
g=np.arange(0,1.025,.025 )
Ls=[500]

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
    

    

fig,ax = plt.subplots(dpi=600)


# Have a look at the colormaps here and decide which one you'd like:
# http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html







Omega = 0.1
gi=[0.2,0.75]
ysq=[squeezednes_i(0.2, L), squeezednes_i(0.75, L)]
fig,ax = plt.subplots(dpi=600)


ax.plot(g, ys01[0], color='lightskyblue', label=r'$\omega_{0}= 0.1 t_{h}$', zorder=1)
ax.plot(g, ys1[0], color='black', label=r'$\omega_{0}= t_{h}$', zorder=1)

ax.plot(g, ys10[0], color='peru', label=r'$\omega_{0}= 10 t_{h}$', zorder = 1)

ax.scatter(gi,ysq,  s=50,marker='x',color='r',linewidths=1.5, zorder=2)


print(double_data(0.2, L), double_data(0.75, L))
    
    



plt.xlabel(r"$g$")
plt.ylabel(r"$\frac{\Delta P}{\Delta X}$")

plt.legend(loc=5, bbox_to_anchor=(0.98, 0.4))



plt.show()



# now make a circle with no fill, which is good for hi-lighting key results





