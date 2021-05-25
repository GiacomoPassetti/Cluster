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
L=500
B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
X= B+Bd
Y= 1j*(Bd-B)
XX=X.dot(X)
YY=Y.dot(Y)
Omega=0.1
g=0.75


dx = double_data(g, L)[0]
dy = double_data(g, L)[1]
print(dx, dy)
Ely = Ellipse((0, 0), dx, dy, color = 'black', fill = False)
circle2 = plt.Circle((0.5, 0.5), 0.4, color='black', fill= False)


fig, ax = plt.subplots(dpi=1000) 
ax.set_aspect('equal')# note we must use plt.subplots, not plt.subplot
# (or if you have an existing figure)
# fig = plt.gcf()
# ax = fig.gca()

ax.arrow(-0.75, 0, 1.5, 0, head_width=0.05, head_length=0.05, fc='k', ec='k')
ax.arrow(0, -0.75, 0, 1.5, head_width=0.05, head_length=0.05, fc='k', ec='k')

ax.add_patch(Ely)
