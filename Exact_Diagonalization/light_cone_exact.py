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
    NS=[]
  
  def expectation(self, Op):
      val=np.real_if_close(np.tensordot(self.v.conj().T, np.tensordot(Op, self.v, 1), 1))
      
      return val      

  def apply(self, Op):
    self.v=np.tensordot(Op, self.v, 1)
  
  def Nf(self):
      nf=[]
      for i in range(L):
          nf.append(self.expectation(N(i+1)))
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
    
def Ground_state_peier(g, J, Omega):
        w, v= eigh(Peier_open(g,Omega, J), eigvals_only=False) 
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
        return  GS
        
    
    

     # This is the function that initializes the fermionic chain in the non-int ground and the Phot
     # in the coherent state, then evolve the system after a small local perturb and save the instantaneous occ.
def light_cone(Omega, J, g, dt, tmax):
    ID='Omega_'+str(Omega)+'J_'+str(J)+' g_'+str(g)+' Nmax_'+str(Nmax)+' L_'+str(L)
    A=Operator_builder([B+Bd]+[Idf]*L)[1]
    pert=Operator_builder([Idb]+[Idf]*int(L/2)+[Cd-0.95*C]+[Idf]*int((L/2)-1))[1]
    U=U_dt(Peier_open(g, Omega, J), dt) #Remember that in order to perform the quench the TE operator 
    GS=Ground_state_peier(0, J, Omega) # and the ground state function must have different g values
    print('Diagonalization done for g='+str(g)+' done')
    GS.apply(displacement(1))
    GS.apply(pert)
    t=list(np.arange(0,tmax, dt))
    n_i_t=[]
    X_t=[]
    for i in t:
        print('begin time step:', i)
        GS.apply(U)
        n_i_t.append(GS.Nf())
        X_t.append(GS.expectation(A))
        print('Done')
        
    np.save('Time_occup_Exact_'+ID, n_i_t)
    np.save('Time_Displ_Exact'+ID, X_t)
    plt.plot(t, X_t)
    plt.xlabel('t')
    plt.ylabel('<X(t)> for g= '+str(g))
    
    
        
    
def Oscillator_Dampening(Omega, J, g, dt, tmax):
    ID='Omega_'+str(Omega)+'J_'+str(J)+' g_'+str(g)+' Nmax_'+str(Nmax)+' L_'+str(L)
    A=Operator_builder([B+Bd]+[Idf]*L)[1]

    U=U_dt(Peier_open(g, Omega, J), dt) #Remember that in order to perform the quench the TE operator 
    GS=Ground_state_peier(0, J, Omega) # and the ground state function must have different g values
    print('Diagonalization done for g='+str(g)+' done')
    GS.apply(displacement(1))

    t=list(np.arange(0,tmax, dt))

    X_t=[]
    for i in t:
        print('begin time step:', i)
        GS.apply(U)

        X_t.append(GS.expectation(A))
        print('Done')
        

    np.save('Time_Displ_Exact_Damp'+ID, X_t)

    return X_t

                      
    
    

Omega, J, g, Nmax, L = 10,1,0,6,8
dt=0.01
tmax=6
filling=int(L/2)
ID='Omega_'+str(Omega)+'J_'+str(J)+' g_'+str(g)+' Nmax_'+str(Nmax)+' L_'+str(L)
Fock=(Nmax+1)*(2**L)
C, Cd, Nf, Idf = FermionSite(None, filling=0.5).C.to_ndarray(), FermionSite(None, filling=0.5).Cd.to_ndarray(), FermionSite(None, filling=0.5).N.to_ndarray(), FermionSite(None, filling=0.5).Id.to_ndarray()
B, Bd, Nb, Idb = BosonSite(Nmax=Nmax,conserve=None, filling=0 ).B.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Bd.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).N.to_ndarray(), BosonSite(Nmax=Nmax,conserve=None, filling=0 ).Id.to_ndarray()
Vac_b=np.zeros((Nmax+1))
Vac_b[0]=1

empty=np.array([1,0])
full=np.array([0,1])
empty_vector=vec_builder([Vac_b]+[empty]*L)[1].reshape(Fock,1)
full_vector=vec_builder([Vac_b]+[full]*L)[1].reshape(Fock,1)
hf=Vector(vec_builder([Vac_b]+[empty]*int(L/2)+ [full]*int(L/2))[1].reshape(Fock, 1))




X=[]
for gg in [0.25, 0.75, 1.5]:
    t=np.arange(0,0.5, 0.1)
    X.append(Oscillator_Dampening(Omega, J, gg, dt, tmax))
t=np.arange(0,tmax,dt)
plt.plot(t,X[0],t,X[1],t,X[2])
plt.xlabel('t')
plt.ylabel(r'$<X(t)>    \Omega=$'+str(Omega))
plt.legend(['g=0.25', 'g=0.75', 'g=1.5'])
    
    