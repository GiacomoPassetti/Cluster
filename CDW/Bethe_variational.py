# -*- coding: utf-8 -*-
"""
Created on Mon May 31 16:15:53 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt
from tenpy.networks.site import BosonSite
from scipy.linalg import expm, sinm, cosm, eigh
from scipy.sparse.linalg import eigsh

Nmax=400
L=1000
B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
X= B+Bd
Y= 1j*(Bd-B)
BB = B.dot(B)
BdBd= Bd.dot(Bd)
GS = np.zeros(Nmax+1)
GS[0] = 1
alpha = 1
eta = 1


def D(alpha):
    D = expm((alpha*Bd)- (np.conj(alpha)*B))
    return D

def S(eta):
    S = expm((1/2)*((np.conj(eta)*BB)-(eta*BdBd)))
    return S

def Ph(g):
    P = expm(1j*(X)*g)
    return P
def Exp(ket, op):
    av = ket.conj().T.dot(op.dot(ket))
    return av


ket = D(alpha).dot(S(eta).dot(GS))

print(abs(Exp(ket, Ph(1))))

alphas = np.arange(0,12.5, 0.5)
etas = np.arange(0, 12.5, 0.5)
a = []
for i in alphas:
    b = []
    for j in etas:
        ket = D(i).dot(S(j).dot(GS))
        b.append(abs(Exp(ket, Ph(0.5))))
    a.append(b)

plt.figure()
plt.imshow(a[::-1],'plasma',
               vmin=None,
               aspect='auto',
               interpolation='nearest',
               extent=(0, 5, 0, 2))

plt.xlabel(r'$U$')
plt.ylabel(r'$g$')


plt.colorbar().set_label(r'$\sum_{i}(n_{i}-\frac{1}{2})^{2}$')




