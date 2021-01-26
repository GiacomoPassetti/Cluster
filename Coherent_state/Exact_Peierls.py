# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 15:28:01 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
from tenpy.networks.site import FermionSite
from tenpy.networks.site import BosonSite


Omega, t, g, Nmax, L = 0,0,0,4,4

def Block(L, op):
    for i in range
def Hamiltonian(Omega, t, g, Nmax, L):
    C, Cd, Nf, Idf = FermionSite(None, filling=0.5).C.to_ndarray(), FermionSite(None, filling=0.5).Cd.to_ndarray(), FermionSite(None, filling=0.5).N.to_ndarray(), FermionSite(None, filling=0.5).Id.to_ndarray()
    B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
    
    for i in range(int(L/2)):
        print(i)
    H=np.zeros((Nmax+2*L+1, Nmax+2*L+1))
 
    return H

A=FermionSite(None, filling=0.5).Cd.to_ndarray()
H=Hamiltonian(Omega, t, g, Nmax, L)
print(A)