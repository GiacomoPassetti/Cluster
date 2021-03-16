# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 15:28:01 2021

@author: giaco
"""

import numpy as np

import scipy.linalg as la
from scipy.linalg import expm, sinm, cosm, eigh
from tenpy.networks.site import FermionSite
from tenpy.networks.site import BosonSite
import copy
import time



  
    
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
    self.v = A
    
    self.norm= np.tensordot(self.v.conj().T, self.v, 1)
    self.NN=Operator_builder([NNb]+[Idf]*L)[1]
    self.NF=[]
    for i in range(L):
        self.NF.append(N(i+1))
     
  
  def expectation(self, Op):
      val=np.real_if_close(np.tensordot(self.v.conj().T, np.tensordot(Op, self.v, 1), 1))
      
      return val      

  def apply(self, Op):
    self.v=np.tensordot(Op, self.v, 1)
  
  def Nf(self):
      nf=[]
      for i in range(L):
          nf.append(self.expectation(self.NF[i]))
      return nf
  def Nb(self):
      return self.expectation(N(0))
  def NNb(self):
      return self.expectation(self.NN)

    
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
    




def peier_periodic(g,Omega,J):
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

    hop_R=[Idb]+[JW]*L
    hop_L=[Idb]+[JW]*L
    hop_R[0]=expm(1j*g*(B+Bd))
    hop_R[L]=-J*Cd
    hop_R[1]=np.tensordot(JW, C, 1)
    hop_L[0]=expm(-1j*g*(B+Bd))
    hop_L[L]=-J*C
    hop_L[1]=np.tensordot(Cd, JW, 1)
    
    kin=kin+ Operator_builder(hop_L)[1]+Operator_builder(hop_R)[1]
    H=kin+cav
    return H

def peier_open(g,Omega,J):
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



def cdw(g,Omega,J, U):
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

    hop_R=[Idb]+[JW]*L
    hop_L=[Idb]+[JW]*L
    hop_R[0]=expm(1j*g*(B+Bd))
    hop_R[L]=-J*Cd
    hop_R[1]=np.tensordot(JW, C, 1)
    hop_L[0]=expm(-1j*g*(B+Bd))
    hop_L[L]=-J*C
    hop_L[1]=np.tensordot(Cd, JW, 1)
    kin=kin+ Operator_builder(hop_L)[1]+Operator_builder(hop_R)[1]
    dN=np.zeros((Fock,Fock))
    for i in range(L-1):
        hop_R=[Idb]+[Idf]*L
        hop_R[i+1]=U*(Nf-(0.5*Idf))
        hop_R[i+2]=(Nf-(0.5*Idf))
        dN=dN+Operator_builder(hop_R)[1]
    hop_R=[Idb]+[Idf]*L
    hop_R[L]=U*(Nf-(0.5*Idf))
    hop_R[1]=(Nf-(0.5*Idf))
    dN=dN+Operator_builder(hop_R)[1]
    
    

    H=kin+cav+dN
    return H      
    




def GS(H):
        w, v= eigh(H, eigvals_only=False) 

        oc=0
        for i in range(Fock):
            print(i)
            GS=Vector(v[:, i])
            energy=w[i]
            oc=sum(GS.Nf())
            if abs(oc-filling)<0.00001:
                break
            else:
                continue

        
        return  GS, energy
    

    
    

Omega, J, g, U, Nmax, L = 10,1,0,0,8,4
h=1
mu=0
dt=0.1
filling=int(L/2)
ID='Omega_'+str(Omega)+'J_'+str(J)+' g_'+str(g)+' Nmax_'+str(Nmax)+' L_'+str(L)+'U_'+str(U)
Fock=(Nmax+1)*(2**L)
C, Cd, Nf, Idf = FermionSite(None, filling=0.5).C.to_ndarray(), FermionSite(None, filling=0.5).Cd.to_ndarray(), FermionSite(None, filling=0.5).N.to_ndarray(), FermionSite(None, filling=0.5).Id.to_ndarray()
B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
NNb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).NN.to_ndarray()
JW=FermionSite(None, filling=0.5).JW.to_ndarray()

"""
ts=time.time()
H=cdw(g, Omega, J, U)
U=expm(-1j*dt*H)
gs, energy=GS(H)
for i in range(10):
    print(time.time()-ts)
    gs.apply(U)
print(energy)
print(gs.Nb())
plt.plot(gs.Nf())
"""
"""
H=cdw(0, Omega, J, 0.2)
gs, energy=GS(H)
print(energy)
"""
# Let's set the Quench protocol:
    #1) Initialize ground for U=0, the ramp will go from 0 to 4
    # we will need bench mark for the energies at some distanced intervals:
"""
ts=time.time()
Us=np.arange(0,4.2,0.2)
gg=np.arange(0,1.2,0.2)
for g in gg:
  for U in Us:
    print('for g'+str(g)+'U'+str(U), time.time()-ts)
    ID='Omega_'+str(Omega)+'J_'+str(J)+' g_'+str(g)+' Nmax_'+str(Nmax)+' L_'+str(L)+'U_'+str(U)
    H=cdw(g, Omega, J, U)
    gs, energy=GS(H)
    np.save(ID+'groundstate.npy', gs.v)
    np.save(ID+'Energy_ground.npy', energy)
    H=0
"""


A=Vector(np.array([1,2,12,1]))



        
       

