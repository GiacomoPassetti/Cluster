# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 15:28:01 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
from scipy.linalg import expm, sinm, cosm
from tenpy.networks.site import FermionSite
from tenpy.networks.site import BosonSite
import copy


Omega, J, g, Nmax, L = 0,1,0,2,2
C, Cd, Nf, Idf = FermionSite(None, filling=0.5).C.to_ndarray(), FermionSite(None, filling=0.5).Cd.to_ndarray(), FermionSite(None, filling=0.5).N.to_ndarray(), FermionSite(None, filling=0.5).Id.to_ndarray()
B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
Vac_b=np.zeros((Nmax+1))
Vac_b[0]=1
empty=np.array([1,0])
full=np.array([0,1])


def Block(A, B):

    C=np.block([
    [A,               np.zeros((A.shape[0], B.shape[0]))],
    [np.zeros((B.shape[0], A.shape[0])), B               ]
    ])
    return C


def Hamiltonian_Peierls(Omega, t, g, Nmax, L):
    C, Cd, Nf, Idf = FermionSite(None, filling=0.5).C.to_ndarray(), FermionSite(None, filling=0.5).Cd.to_ndarray(), FermionSite(None, filling=0.5).N.to_ndarray(), FermionSite(None, filling=0.5).Id.to_ndarray()
    B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
    
    for i in range(int(L/2)):
        print(i)
    H=np.zeros((Nmax+2*L+1, Nmax+2*L+1))
 
    return H


def Operator_builder(List):
    Op=np.tensordot(List[0], List[1], axes=0)
    ind=list(np.arange(0,2*(L+1),2))+list(np.arange(1,2*(L+1),2))
    
    for i in range(L-1):
       Op=np.tensordot(Op, List[i+2], axes=0)
    Op=np.transpose(Op, ind)
    Opq=np.reshape(Op, [(Nmax+1)*(2**L), (Nmax+1)*(2**L)])
    return Op, Opq

"""
H=np.zeros(((Nmax+1)*(2**L), (Nmax+1)*(2**L)))
for i in range(L-1):
     I=[Idb]+[Idf]*L
     I[0]=expm(1j*g*(B+Bd))
     I[i+1]=-J*Cd
     I[i+2]=C
     H=H+Operator_builder(I)
H=H+H.conj().T  
"""
 
     
        
class Vector:
  def __init__(self, A):
    self.vector = A
     
  
  def expectation(self, Op):
      val=np.tensordot(self.vector.conj().T, np.tensord(Op, self.vector, 1), 1)
      return val      

  def apply(self, Op):
    self.vector=np.tensordot(Op, self.vector, 1)
    
  def Vector_builder(self, List):
    v=np.tensordot(List[0], List[1], axes=0)
    
    
    for i in range(L-1):
       v=np.tensordot(v, List[i+2], axes=0)
    print(v.shape)
    v=np.reshape(v, [(Nmax+1)*(2**L)])
    print(v.shape)
    self.vector=v
     



       
       
       
       

