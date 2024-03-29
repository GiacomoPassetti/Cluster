# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 20:09:40 2021

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
from tenpy.linalg.charges import LegCharge, ChargeInfo
from tenpy.algorithms.truncation import truncate, svd_theta, TruncationError
import tenpy.linalg.np_conserved as npc
from scipy.linalg import expm
import pickle
import matplotlib.pyplot as plt
import time

def sites(L,Nmax):
 FSite=FermionSite('N', filling=0.5)
 qflat=[[0]]*(Nmax+1)
 ch=ChargeInfo([1], names=None)
 leg = LegCharge.from_qflat(ch, qflat, qconj=1)
 BSite=BosonSite(Nmax=Nmax,conserve=None, filling=0 )
 BSite.change_charge(leg)
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

def last_site(L):
    ps=['vac']
    for i in range(int(L-1)):
        ps.append('empty')
    ps.append('full')
    return ps

def first_site(L):
    ps=['vac']
    ps.append('full')
    for i in range(int(L-2)):
        ps.append('empty')
    
    ps.append('full')
    return ps
def left(L):
    ps=['vac']
    
    for i in range(int(L/2)):
        ps.append('full')
    for i in range(int(L/2)):
        ps.append('empty')

    return ps
def mixed_state(L):
    ps=['vac']
    ms = np.array([1/np.sqrt(2), 1/np.sqrt(2)])
    for i in range(int(L/2)):
        ps.append(ms)
        ps.append(ms)
    return ps

def H_Peier_bond(psi, g, J, Omega, V, h1, h2, L):
    
    #In order to read quickly the total energy I define both the bond energy for the coupling with the cavity and for only the fermions
    Peier=npc.outer(npc.expm(1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']),npc.outer(-J*psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].C.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
    Peier_hc=npc.outer(npc.expm(-1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']), npc.outer(-J*psi.sites[1].C.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
    cav=npc.outer((Omega/((L-1)))*psi.sites[0].N.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
    ons_l=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(h1*psi.sites[1].N.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
    rep=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(V*psi.sites[1].N.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].N.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
    ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),h2*psi.sites[1].N.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
    H_bond=Peier+Peier_hc+cav+ons_l+ons_r+rep  #This is the energetic term that will be used in the TEBD algorithm
    return  H_bond

def H_Peier_bond_dn(psi, g, J, Omega, V, L):
    
    #In order to read quickly the total energy I define both the bond energy for the coupling with the cavity and for only the fermions
    Peier=npc.outer(npc.expm(1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']),npc.outer(-J*psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].C.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
    Peier_hc=npc.outer(npc.expm(-1j*g*(psi.sites[0].B+psi.sites[0].Bd)).replace_labels(['p', 'p*'], ['p0', 'p0*']), npc.outer(-J*psi.sites[1].C.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Cd.replace_labels(['p', 'p*'], ['p2', 'p2*']))).itranspose([0,2,4,1,3,5])
    cav=npc.outer((Omega/((L-1)))*psi.sites[0].N.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
    
    rep=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(V*psi.sites[1].dN.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].dN.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
    
    H_bond=Peier+Peier_hc+cav+rep  #This is the energetic term that will be used in the TEBD algorithm
    return  H_bond


def U_bond(dt, H_bond):
    
    H2 = H_bond.combine_legs([('p0', 'p1', 'p2'), ('p0*', 'p1*', 'p2*')], qconj=[+1, -1])
    H2= (-dt)*H2
    U=npc.expm(H2).split_legs()
    return U

def from_full_custom(
                  siti,
                  theta,
                  trunc_par,
                  outer_S,
                  cutoff=1.e-16,
                  form=None,
                  normalize=True,
                  ):

        
        L = len(siti)

        B_list = [None] * L
        S_list = [None] * (L + 1)
        norm = 1. 
        
        labels = ['vL'] + ['p' + str(i) for i in range(L)] + ['vR']
        theta.itranspose(labels)
        # combine legs from left
        for i in range(0, L - 1):
            theta = theta.combine_legs([0, 1])  # combines the legs until `i`
        # now psi has only three legs: ``'(((vL.p0).p1)...p{L-2})', 'p{L-1}', 'vR'``
        for i in range(L - 1, 0, -1):
            # split off B[i]
            theta = theta.combine_legs([labels[i + 1], 'vR'])
            theta, S, B, err, renorm = svd_theta(theta, trunc_par, qtotal_LR=[None, None], inner_labels=['vR', 'vL'])
            
            
            if i > 1:
                theta.iscale_axis(S, 1)
            B_list[i] = B.split_legs(1).replace_label(labels[i + 1], 'p')
            S_list[i] = S
            theta = theta.split_legs(0)
        # psi is now the first `B` in 'A' form
        B_list[0] = theta.replace_label(labels[1], 'p')
        B_form = ['A'] + ['B'] * (L - 1)
        S_list[0], S_list[-1] = outer_S
        res = MPS(siti, B_list, S_list, bc='segment', form=B_form, norm=norm)
        if form is not None:
            res.convert_form(form)
        return res, err




def ansatz_wf(Nmax, L):
    ps= product_state(L)
    site= sites(L,Nmax)
    psi=MPS.from_product_state(site, ps)
    return psi

def ansatz_last(Nmax, L):
    ps=last_site(L)
    site= sites(L,Nmax)
    psi=MPS.from_product_state(site, ps)
    return psi

def ansatz_first(Nmax, L):
    ps=first_site(L)
    site= sites(L,Nmax)
    psi=MPS.from_product_state(site, ps)
    return psi

def ansatz_left(Nmax, L):
    ps=left(L)
    site= sites(L,Nmax)
    psi=MPS.from_product_state(site, ps)
    return psi

def apply_local_cav_r(psi, i, op, trunc_param):
            
            cutoff=1.e-13
            "1--  Applico U"

            
            # Prendo il tensore di rango 3 che mi serve e contraggo
            n = 3
            p = psi._get_p_labels(n, False)
            pstar = psi._get_p_labels(n, True)
            th = psi.get_theta(i, n)
            th = npc.tensordot(op, th, axes=[pstar, p])
            
            "2-- Permutazione di indici che realizza lo swap"
            
            th.ireplace_labels(['p0', 'p1','p2'], ['p1', 'p0', 'p2'])
            
            
            
            "3-- Scomposizione in A S B B e ridefinizione di psi"
            split_th, err = from_full_custom(psi.sites[i:i + n], th, trunc_param,outer_S= (psi.get_SL(i), psi.get_SR(i + n - 1)))
            for j in range(n):
                psi.set_B(i + j, split_th._B[j], split_th.form[j])
            for j in range(n - 1):
                psi.set_SR(i + j, split_th._S[j + 1])
            siteL, siteR = psi.sites[psi._to_valid_index(i)], psi.sites[psi._to_valid_index(i + 1)]
            psi.sites[psi._to_valid_index(i)] = siteR  # swap 'sites' as well
            psi.sites[psi._to_valid_index(i + 1)] = siteL
            
            return err
            


def apply_local_cav_end(psi, i, op, trunc_param):
            cutoff=1.e-13
            "1--  Applico U"

            
            # Prendo il tensore di rango 3 che mi serve e contraggo
            n = 3
            p = psi._get_p_labels(n, False)
            pstar = psi._get_p_labels(n, True)
            th = psi.get_theta(i, n)
            th = npc.tensordot(op, th, axes=[pstar, p])
            

            "3-- Scomposizione in A S B B e ridefinizione di psi"
            split_th, err = from_full_custom(psi.sites[i:i + n], th, trunc_param,outer_S= (psi.get_SL(i), psi.get_SR(i + n - 1)))
            for j in range(n):
                psi.set_B(i + j, split_th._B[j], split_th.form[j])
            for j in range(n - 1):
                psi.set_SR(i + j, split_th._S[j + 1])
            
            return err


def apply_local_cav_l(psi, i, op, trunc_param):
            i=i-1
            cutoff=1.e-13
            "1--  Genero theta e swappo left"
            n=3
            th = psi.get_theta(i, n)
            th.ireplace_labels(['p0', 'p1','p2'], ['p1', 'p0', 'p2'])
            
            p = psi._get_p_labels(n, False)
            pstar = psi._get_p_labels(n, True)

            
            "2-- applico"

            
            th = npc.tensordot(op, th, axes=[pstar, p])
            
            
            
            
            "3-- Scomposizione in A S B B e ridefinizione di psi"
            split_th, err = from_full_custom(psi.sites[i:i + n], th, trunc_param,outer_S= (psi.get_SL(i), psi.get_SR(i + n - 1)))
            for j in range(n):
                psi.set_B(i + j, split_th._B[j], split_th.form[j])
            for j in range(n - 1):
                psi.set_SR(i + j, split_th._S[j + 1])
            siteL, siteR = psi.sites[psi._to_valid_index(i)], psi.sites[psi._to_valid_index(i + 1)]
            psi.sites[psi._to_valid_index(i)] = siteR  # swap 'sites' as well
            psi.sites[psi._to_valid_index(i + 1)] = siteL
            
            return err
            

def Energy(psi, H_bond, L, trunc_param):
         
        E=[]

        for i in range(int(L/2)-1): # First Odd sweep
            
            psi.swap_sites(2*i, swap_op=None, trunc_par=trunc_param)
            E.append(psi.expectation_value(H_bond[2*i+1], [2*i+1]))
            psi.swap_sites(2*i+1, swap_op=None, trunc_par=trunc_param)

        for i in range(int(L/2)-1):
            

            E.append(psi.expectation_value(H_bond[-2*i-1], [L-2-2*i]))
            psi.swap_sites(L-3-2*i, swap_op=None, trunc_par=trunc_param)
            psi.swap_sites(L-4-2*i, swap_op=None, trunc_par=trunc_param)

        E.append(psi.expectation_value(H_bond[0], [0]))
        E_tot=np.sum(E)

        
        return E_tot




def full_sweep(psi, step, U, Id, trunc_param, L):
  eps=TruncationError()

  for _ in range(step):
      

   for i in range((L-2)//2):
     
     eps += apply_local_cav_r(psi, 2*i, U[2*i], trunc_param) 
     eps += apply_local_cav_r(psi, 2*i+1, Id, trunc_param)  
   
   eps += apply_local_cav_end(psi, L-2, U[L-2], trunc_param)
   for i in range((L-2)//2):
    
    eps += apply_local_cav_l(psi, L-2-2*i, U[-2-2*i], trunc_param)
    eps += apply_local_cav_l(psi, L-3-2*i, Id, trunc_param)
  return eps

    


def Suz_trot_im(psi, delta_t, max_error_E, N_steps, H_bond, trunc_param, L, Id):
 start_time=time.time()
 DeltaE=2*max_error_E
 E_old=Energy(psi, H_bond, L, trunc_param)
 for dt in range(len(delta_t)):
    print("delta_tau =", delta_t[dt], "Time of evaluation:", time.time()-start_time)

    U=[]
    for i in range(L-1):
        U.append(U_bond(delta_t[dt], H_bond[i]))
    DeltaE= 2*max_error_E[dt]
    step=0
    while (DeltaE > max_error_E[dt]):
    
      
      full_sweep(psi, N_steps[dt], U, Id, trunc_param, L)
      print(psi.correlation_function('N', 'N', sites1=[1], sites2=range(1, L+1)))
      plt.plot(psi.expectation_value('N'))
      plt.show()
        

      step += N_steps[dt]
      E=Energy(psi, H_bond, L, trunc_param)
      DeltaE=np.abs(E_old-E)
      E_old=E
      if step > 1200:
          break

      
      print("After", step, "steps, E_tot = ", E, "and DeltaE = ", DeltaE , "Time of evaluation:", time.time()-start_time)



def Iterative_g(psi_f, g1, g2, step, L, Omega, J, h, V, Nmax):
    print(psi_f)
    siti = sites(L,Nmax)
    ps= product_state(L)
    N_steps=[10, 10, 10, 10, 10, 10]
    delta_t_im=[0.1, 1.e-2, 1.e-3, 1.e-4, 1.e-5]
    trunc_param={'chi_max':120,'svd_min': 1.e-13, 'verbose': False}
    max_error_E=[1.e-8, 1.e-7, 1.e-6, 1.e-6, 1.e-6, 1.e-6]
    gs=np.arange(g1,g2+step, step)
    
    # Produce the first psi for g very small
    
    ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(0.1)
    g=0.5/np.sqrt(L)
    psi=MPS.from_product_state(siti, ps)
    for i in range(L):
          psi.set_B(i+1, psi_f.get_B(i))
          psi.set_SL(i+1, psi_f.get_SL(i))
          psi.set_SR(i+1, psi_f.get_SR(i))
          
          
    
    Id=ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
    H_bond=[]
    for i in range(L-1):
          H_bond.append(H_Peier_bond(psi, g, J, Omega,V, (2*i+1)*h, (2*i+2)*h, L))
   
    Suz_trot_im(psi, delta_t_im, max_error_E, N_steps, H_bond, trunc_param, L, Id)
    with open(ID+'imTEBD.pkl', 'wb') as f:
          pickle.dump(psi, f)
          
          # Then iterate starting from the Ansatz generated
          
    N_steps=[ 10]
    delta_t_im=[ 1.e-5]
    trunc_param={'chi_max':120,'svd_min': 1.e-13, 'verbose': False}
    max_error_E=[ 1.e-6]
    for g0 in gs:
        ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g0)
        g=g0/np.sqrt(L)
        H_bond=[]
        for i in range(L-1):
          H_bond.append(H_Peier_bond(psi, g, J, Omega,V, (2*i+1)*h, (2*i+2)*h, L))
        Suz_trot_im(psi, delta_t_im, max_error_E, N_steps, H_bond, trunc_param, L, Id)
        with open(ID+'imTEBD.pkl', 'wb') as f:
          pickle.dump(psi, f)


def Iterative_g_from_load(psi, g1, g2, step, L, Omega, J, h, V, Nmax, N_steps, delta_t_im, trunc_param, max_error):
    
    siti = sites(L,Nmax)
    ps= product_state(L)
    N_steps=[10, 10, 10, 10, 10, 10]
    delta_t_im=[0.1, 1.e-2, 1.e-3, 1.e-4, 1.e-5]
    trunc_param={'chi_max':120,'svd_min': 1.e-13, 'verbose': False}
    max_error_E=[1.e-8, 1.e-7, 1.e-6, 1.e-6, 1.e-6, 1.e-6]
    gs=np.arange(g1,g2+step, step)
    
    # Produce the first psi for g very small
    
    ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(0.1)
    g=0.1/np.sqrt(L)
    
    
    Id=ons_r=npc.outer(psi.sites[0].Id.replace_labels(['p','p*'],['p0', 'p0*']),npc.outer(psi.sites[1].Id.replace_labels(['p', 'p*'], ['p1', 'p1*']),psi.sites[1].Id.replace_labels(['p', 'p*'], ['p2', 'p2*'])) ).itranspose([0,2,4,1,3,5])
    """
    H_bond=[]
    
    for i in range(L-1):
          H_bond.append(H_Peier_bond(psi, g, J, Omega,V, (2*i+1)*h, (2*i+2)*h, L))
   
    Suz_trot_im(psi, delta_t_im, max_error_E, N_steps, H_bond, trunc_param, L, Id)
    with open(ID+'imTEBD.pkl', 'wb') as f:
          pickle.dump(psi, f)
          
          # Then iterate starting from the Ansatz generated
    """
    N_steps=[ 10]
    delta_t_im=[ 1.e-5]
    trunc_param={'chi_max':120,'svd_min': 1.e-13, 'verbose': False}
    max_error_E=[ 1.e-6]
    for g0 in gs:
        ID='Psi_GS_Nmax_'+str(Nmax)+'L_'+str(L)+'Omega_'+str(Omega)+'J_'+str(J)+'h_'+str(h)+'V_'+str(V)+'g_0'+str(g0)
        g=g0/np.sqrt(L)
        H_bond=[]
        for i in range(L-1):
          H_bond.append(H_Peier_bond(psi, g, J, Omega,V, (2*i+1)*h, (2*i+2)*h, L))
        Suz_trot_im(psi, delta_t_im, max_error_E, N_steps, H_bond, trunc_param, L, Id)
        with open(ID+'imTEBD.pkl', 'wb') as f:
          pickle.dump(psi, f)


