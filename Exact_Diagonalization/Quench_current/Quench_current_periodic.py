# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 10:02:17 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
from scipy.linalg import expm, sinm, cosm, eigh
from tenpy.networks.site import FermionSite
from tenpy.networks.site import BosonSite
import copy



  
    
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
    
def displacement(alpha):
    l=[expm(alpha*Bd-np.conj(alpha)*B)]+[Idf]*L
    l=Operator_builder(l)[1]
    return l
    
def U_dt(op, dt):
    U=expm(-1j*dt*op)
    return U


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

def Peier_impurity(g,Omega,J, eta):
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
    impurity=[Idb]+[Idf]*L
    impurity[1]=eta*Nf
    imp=Operator_builder(impurity)[1]
    H=kin+cav+imp
    return H

def Peier_current(g,Omega,J, dV):
    cav=Operator_builder([Omega*Nb]+[Idf]*L)[1]
    kin=np.zeros((Fock,Fock))
    ons=np.zeros((Fock,Fock))
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
    for i in range(L):
        II=[Idb]+[Idf]*L
        II[i+1]=(-dV/(L-1)*i+(dV/2))*Nf
        ons=ons+Operator_builder(II)[1]
    
    H=kin+cav+ons
    return H

def Peier_periodic(g,Omega,J):
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
    hop_R=[expm(1j*g*(B+Bd))]+[C]+[Idf]*(L-2)+[-J*Cd]
    hop_L=[expm(-1j*g*(B+Bd))]+[-J*Cd]+[Idf]*(L-2)+[C]
    kin=kin+ Operator_builder(hop_L)[1]+Operator_builder(hop_R)[1]
    H=kin+cav
    return H

def Peier_current_periodic(g,Omega,J, dV):
    cav=Operator_builder([Omega*Nb]+[Idf]*L)[1]
    kin=np.zeros((Fock,Fock))
    ons=np.zeros((Fock,Fock))
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
    hop_R=[expm(1j*g*(B+Bd))]+[C]+[Idf]*(L-2)+[-J*Cd]
    hop_L=[expm(-1j*g*(B+Bd))]+[-J*Cd]+[Idf]*(L-2)+[C]
    kin=kin+ Operator_builder(hop_L)[1]+Operator_builder(hop_R)[1]
    for i in range(L):
        II=[Idb]+[Idf]*L
        II[i+1]=(-dV/(L-1)*i+(dV/2))*Nf
        ons=ons+Operator_builder(II)[1]
    
    H=kin+cav+ons
    return H
    
    
def plot_ph_av(gmin, gmax, steps, Omega, J):
    fot_avg=[]
    err=filling
    gs=list(np.arange(gmin, gmax, steps))
    
    for g in gs:
        
        w, v= eigh(Peier_open(g,Omega, J), eigvals_only=False) 
        vectors=[]
        oc=0
        for i in range(Fock):
            print(i)
            vectors.append(Vector(v[:, i]))
            oc=sum(vectors[i].Nf())
            GS=vectors[i]
            E_gs=w[i]
            if abs(oc-filling)<0.00001:
                break
            else:
                continue

        fot_avg.append(GS.Nb())
    plt.plot(gs, fot_avg)
    plt.xlabel('g')
    plt.ylabel(r'$<N_{ph}>$')
    plt.show
    np.save('C:/users/giaco/Desktop/Cluster/Exact_Diagonalization/Data/N_avg_bosEXACT'+ID, fot_avg)
    
def Ground_state_peier_periodic(g, J, Omega):
        w, v= eigh(Peier_periodic(g,Omega, J), eigvals_only=False) 
        vectors=[]
        oc=0
        for i in range(Fock):
            
            vectors.append(Vector(v[:, i]))
            oc=sum(vectors[i].Nf())
            GS=vectors[i]
            E_gs=w[i]
            if abs(oc-filling)<0.00001:
                break
            else:
                continue
        ni=GS.Nf()
        plt.plot(ni)
        plt.show()
        
        return  GS
        
    

    

Omega, J, g, Nmax, L = 5,1,1,6,8
eta=0.2
JW=FermionSite(None, filling=0.5).JW.to_ndarray()
filling=int(L/2)
ID='Omega_'+str(Omega)+'J_'+str(J)+' g_'+str(g)+' Nmax_'+str(Nmax)+' L_'+str(L)
Fock=(Nmax+1)*(2**L)
C, Cd, Nf, Idf = FermionSite(None, filling=0.5).C.to_ndarray(), FermionSite(None, filling=0.5).Cd.to_ndarray(), FermionSite(None, filling=0.5).N.to_ndarray(), FermionSite(None, filling=0.5).Id.to_ndarray()
B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
Vac_b=np.zeros((Nmax+1))
Vac_b[0]=1
dV=5
empty=np.array([1,0])
full=np.array([0,1])
empty_vector=vec_builder([Vac_b]+[empty]*L)[1].reshape(Fock,1)
full_vector=vec_builder([Vac_b]+[full]*L)[1].reshape(Fock,1)
hf=Vector(vec_builder([Vac_b]+[empty]*int(L/2)+ [full]*int(L/2))[1].reshape(Fock, 1))






def Quench_current(Omega, J, g, Nmax, L, dV): # Let's see what happens to the cavity when i try to put a current on the wire
    ID='Omega_'+str(Omega)+'J_'+str(J)+' g_'+str(g)+' Nmax_'+str(Nmax)+' L_'+str(L)+'eta_'+str(eta)+ 'dV_'+str(dV)
    GS=Ground_state_peier_periodic(g, J, Omega)
    A=Operator_builder([B+Bd]+[Idf]*L)[1]
    At=[GS.expectation(A)]
    Uimp=U_dt(Peier_periodic(g, Omega, J), 0.05)
    Upert=U_dt(Peier_current_periodic(g,Omega,J, dV), 0.05)
    t1=np.arange(0.05, 1.05, 0.05)
    t2=np.arange(1.05, 50.05, 0.05)
    t=np.arange(0,50.05, 0.05)
    print(len(t1), len(t2))
    for i in list(t1):
       print('Time step:', i)
       GS.apply(Uimp)
       At.append(GS.expectation(A))
    for i in list(t2):
       print('Time step:', i)
       GS.apply(Upert)
       At.append(GS.expectation(A))
    
    np.save('Efield_quench_current'+ID, At)
    plt.plot(t,At)
    plt.xlabel('t')
    plt.ylabel('A(t)')
    plt.show()

Quench_current(Omega, J, g, Nmax, L, dV)
        