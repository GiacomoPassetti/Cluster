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








def H_Peier(J, g, Omega):
    H_peier= Operator_builder([Omega*Nb]+[Idf]*L)[1]
    
    
    
    
    for i in range(L-2):
        hop=[expm(1j*g(B+Bd))]+[Idf]*L
        hop[i+1]=J*Cd
        hop[i+2]=C
        H_peier= H_peier + Operator_builder(hop)[1]
        print(Operator_builder(hop)[1])
    H_peier=H_peier+H_peier.conj().T
    return H_peier

    
  
    
def vec_builder(List):
    vec=np.tensordot(List[0], List[1], axes=0)
    for i in range(L-1):
       vec=np.tensordot(vec, List[i+2], axes=0)
    vecq=np.reshape(vec, ([Fock]))
    return vec, vecq
    

def Operator_builder(List):  # Realizes the tensor dot of the List of operators in entry and can return both the tensorized and the squared form
    Op=np.tensordot(List[0], List[1], axes=0)
    ind=list(np.arange(0,2*(L+1),2))+list(np.arange(1,2*(L+1),2))
    
    for i in range(L-1):
       Op=np.tensordot(Op, List[i+2], axes=0)
    Op=np.transpose(Op, ind)
    Opq=np.reshape(Op, [(Nmax+1)*(2**L), (Nmax+1)*(2**L)])
    return Op, Opq


    
        
class Vector:
  def __init__(self, A):
    self.vector = A
     
  
  def expectation(self, Op):
      val=np.tensordot(self.vector.conj().T, np.tensordot(Op, self.vector, 1), 1)
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
     

Omega, J, g, Nmax, L = 0,1,0,8,8
Fock=(Nmax+1)*(2**L)
C, Cd, Nf, Idf = FermionSite(None, filling=0.5).C.to_ndarray(), FermionSite(None, filling=0.5).Cd.to_ndarray(), FermionSite(None, filling=0.5).N.to_ndarray(), FermionSite(None, filling=0.5).Id.to_ndarray()
B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
Vac_b=np.zeros((Nmax+1))
Vac_b[0]=1
empty=np.array([1,0])
full=np.array([0,1])
empty_vector=vec_builder([Vac_b]+[empty]*L)[1].reshape(Fock,1)
full_vector=vec_builder([Vac_b]+[full]*L)[1].reshape(Fock,1)
half_filled=Vector(vec_builder([Vac_b]+[empty]*int(L/2)+ [full]*int(L/2))[1].reshape(Fock, 1))

def N(i):
    if i==0:
        ni=Operator_builder([Nb]+[Idf]*L)[1]
    else :
        ni=[Idb]+[Idf]*L
        ni[i]=Nf
        ni=Operator_builder(ni)[1]
    return ni

def Nf_tot():
    Nt=np.zeros((Fock,Fock))
    for i in range(L):
        Nt=Nt+N(i+1)
    return Nt
    


def Peier_open(g,Omega,J):
    cav=Operator_builder([Omega*Nb]+[Idf]*L)[1]
    kin=np.zeros((Fock,Fock))
    for i in range(L-1):
        hop_R=[Idb]+[Idf]*L
        hop_L=[Idb]+[Idf]*L
        
        hop_R[0]=expm(1j*g*(B+Bd))
        hop_R[i+1]=-J*Cd
        hop_R[i+2]=C
        hop_L[0]=expm(-1j*g*(B+Bd))
        hop_L[i+1]=-J*C
        hop_L[i+2]=Cd
        kin=kin+ Operator_builder(hop_L)[1]+Operator_builder(hop_R)[1]
    H=kin+cav
    return H

w, v= eigh(Peier_open(g,Omega, J), eigvals_only=False) 





        
       

