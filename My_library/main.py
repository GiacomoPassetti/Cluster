# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 14:42:20 2021

@author: giaco
"""
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/My_library')


import tenpy
import copy
import sys
import numpy as np
import numpy.linalg as alg
from tenpy import models
from tenpy.networks.site import FermionSite
from tenpy.networks.site import BosonSite
from tenpy.networks.mps import MPS
from tenpy.tools.params import get_parameter
import tenpy.linalg.np_conserved as npc
from scipy.linalg import expm
import pickle
import matplotlib.pyplot as plt
import time



def sites(L,Nmax):
 FSite=FermionSite(None, filling=0.5)
 BSite=BosonSite(Nmax=Nmax,conserve=None, filling=0 )
 sites=[]
 sites.append(BSite)
 for i in range(L):
     sites.append(FSite)
 return sites

def product_state(L):
    ps=['vac']
    for i in range(int(L/2)):
        ps.append('empty')
        ps.append('full')
    return ps

def mixed_state(L):
    ps=['vac']
    ms = np.array([1/np.sqrt(2), 1/np.sqrt(2)])
    for i in range(int(L/2)):
        ps.append(ms)
        ps.append(ms)
    return ps


    
def U_bond(dt, H_bond):
    
    H2 = H_bond.combine_legs([('p0', 'p1', 'p2'), ('p0*', 'p1*', 'p2*')], qconj=[+1, -1])
    H2= (-dt)*H2
    U=npc.expm(H2).split_legs()
    return U

def Energy(psi, H_bond, L):
         
        E=[]

        for i in range(int(L/2)-1): # First Odd sweep
            
            psi.swap_sites(2*i, swap_op=None, trunc_par=trunc_param[1])
            E.append(psi.expectation_value(H_bond[2*i+1], [2*i+1]))
            psi.swap_sites(2*i+1, swap_op=None, trunc_par=trunc_param[1])

        for i in range(int(L/2)-1):
            

            E.append(psi.expectation_value(H_bond[-2*i-1], [L-2-2*i]))
            psi.swap_sites(L-3-2*i, swap_op=None, trunc_par=trunc_param[1])
            psi.swap_sites(L-4-2*i, swap_op=None, trunc_par=trunc_param[1])

        E.append(psi.expectation_value(H_bond[0], [0]))
        E_tot=np.sum(E)

        
        return E_tot


def H_Peier_bond(g, J, Omega, h1, h2, L):
    
    #In order to read quickly the total energy I define both the bond energy for the coupling with the cavity and for only the fermions
    Peier=npc.outer(npc.expm(1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']),npc.outer(-J*psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].C.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
    Peier_hc=npc.outer(npc.expm(-1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']), npc.outer(-J*psi.sites[1].C.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
    cav=npc.outer((Omega/((L-1)))*psi.sites[0].N.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
    ons_l=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(h1*psi.sites[1].N.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
    ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),h2*psi.sites[1].N.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
    H_bond=Peier+Peier_hc+cav+ons_l+ons_r  #This is the energetic term that will be used in the TEBD algorithm
    return  H_bond

def Suz_trot_im(psi, delta_t, max_error_E, N_steps, H_bond):
 start_time=time.time()
 DeltaE=2*max_error_E
 E_old=Energy(psi, H_bond, L)
 for dt in range(len(delta_t)):
    print("delta_tau =", delta_t[dt], "Time of evaluation:", time.time()-start_time)

    U=[]
    for i in range(L-1):
        U.append(U_bond(delta_t[dt], H_bond[i]))
    DeltaE= 2*max_error_E[dt]
    step=0
    while (DeltaE > max_error_E[dt]):
    
      for T in range(N_steps[dt]): 
        print("step :", T, "Time of evaluation:", time.time()-start_time)
        for i in range(int(L/2)-1): # First Odd sweep

      
            psi.swap_sites(2*i, swap_op=None, trunc_par=trunc_param[dt])


            psi.apply_local_op((2*i)+1 , U[2*i+1], unitary=False, renormalize=True)


            psi.swap_sites(2*i+1, swap_op=None, trunc_par=trunc_param[dt])
        
        for i in range(int(L/2)-1):# Even sweep
            

            psi.apply_local_op(L-2-2*i, U[-1-2*i], unitary=False, renormalize=True)
            psi.swap_sites(L-3-2*i, swap_op=None, trunc_par=trunc_param[dt])
            psi.swap_sites(L-4-2*i, swap_op=None, trunc_par=trunc_param[dt])


        psi.apply_local_op(0, U[0], unitary=False, renormalize=True)


        psi.compress_svd(trunc_param[dt])
        print(sum(psi.expectation_value('N')))
        

      step += N_steps[dt]
      E=Energy(psi, H_bond, L)
      DeltaE=np.abs(E_old-E)
      E_old=E
      plt.plot(psi.expectation_value('N', [1,2,3,4,5,6]))
      plt.show()
      
      print("After", step, "steps, E_tot = ", E, "and DeltaE = ", DeltaE )




Nmax=8
L=6
g= 0
Omega  = 0
J=1
h=1
V=0
max_error_E=[0.00000001, 1.e-6, 1.e-6, 1.e-7, 1.e-8]
ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)
N_steps=[10, 10, 10,10, 10]
sites = sites(L,Nmax)

ps= product_state(L)
psi=MPS.from_product_state(sites, ps)
cav=npc.outer((Omega/((L-1)))*psi.sites[0].N.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
H_bond=[]
for i in range(L-1):
    H_bond.append(H_Peier_bond(g, J, Omega, (1+i)*h, (2+i)*h, L))

ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'g_'+str(g)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)
delta_t_im=[0.1, 0.1, 0.1, 1.e-4, 1.e-5]
chis=[80, 80, 80, 80, 80]
verbose=False
trunc_param=[]
for i in range(len(chis)):
    trunc_param.append({'chi_max':chis[i],'svd_min': 1.e-13, 'verbose': verbose})
    


Suz_trot_im(psi, delta_t_im, max_error_E, N_steps, H_bond)
plt.plot(psi.expectation_value('N',[1,2,3,4,5,6]))

