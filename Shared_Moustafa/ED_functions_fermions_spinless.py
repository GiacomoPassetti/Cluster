#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 19:17:25 2021

@author: moustafa
"""
import numpy as np
import itertools
import math
from scipy import sparse
from scipy.sparse.linalg import expm
from scipy.sparse.linalg import eigsh 
import matplotlib.pyplot as plt 
from scipy.linalg import cosm


# OBC  = 0 
# PBC  = 1
def states_gen(L,N):
    which = np.array(list(itertools.combinations(range(L), N)))
    #print(which)
    grid = np.zeros((len(which), L), dtype="int8")

    # Magic
    grid[np.arange(len(which))[None].T, which] = 1
    
    return grid

def SzSz_L_red(BC, L, N, N_ph, OBC, PBC):
    dim      = np.int(math.factorial(L)/math.factorial(L-N)/math.factorial(N)*(N_ph+1))
    Id_M_L   = sparse.identity(dim)
    Mrx   = (dim, dim)
    init=sparse.csr_matrix(Mrx)
    for i in range(L-1):
        init += (c_dag_c_i_red(L, N, N_ph, i)-(0.5*Id_M_L)).dot(c_dag_c_i_red(L, N, N_ph, i+1)-(0.5*Id_M_L))
    init += (c_dag_c_i_red(L, N, N_ph, 0)-(0.5*Id_M_L)).dot(c_dag_c_i_red(L, N, N_ph, L-1)-(0.5*Id_M_L))

    return init
    
def kin_L_red (BC, L, N, N_ph, g, OBC, PBC):
    
    states = states_gen(L , N) 
    #print(states)
    num_rows, num_cols = states.shape
    #print(states.shape)
    c_dag_c_kin_L = np.zeros((num_rows , num_rows))

    states_new = states.copy()
    #print(states_new)

    for i in range(num_rows): 
        #print('i is', i)
        for j in range(L - 1):
            #print('j is', j)
            if states_new[i][j] == 0 and states_new[i][j + 1] == 1:
                #print('old',states_new[i][:])
                states_new[i][j]     = 1
                states_new[i][j + 1] = 0
                #print('new', states_new[i][:])

                

                for k in range(num_rows):
                    
                    if np.array_equal(states_new[i][:],states[k][:]):
                        #print('that state is', states_new[i][:])
                        #print('coordinates are', i , k)
                        c_dag_c_kin_L[i][k] = 1
                        break
                
                states_new[i][:] = states[i][:]
                #print('new to old', states_new[i][:])
                #print('old', states[i][:])
        
        if BC == PBC:
            
            if states_new[i][L - 1] == 0 and states_new[i][0] == 1:
                
                states_new[i][0]         = 0  
                
                factor = np.sum(states_new[i][:])
                factor = (-1)**factor 
                
                states_new[i][L - 1]     = 1

                
                for k in range(num_rows):
                    if np.array_equal(states_new[i][:],states[k][:]):
                        #print('that state is', states_new[i][:])
                        #print('coordinates are', i , k)
                        c_dag_c_kin_L[i][k] = factor
                        break
                
                states_new[i][:] = states[i][:]
                #print('new to old', states_new[i][:])
                #print('old', states[i][:])
                
    
    Mrx_a_dag   = Mrces_a(N_ph)[0]
    Mrx_a       = Mrces_a(N_ph)[1]
    
    a_dag_plus_a = Mrx_a_dag + Mrx_a
    
    exp_iA = expm(1j * g * a_dag_plus_a)
    #exp_iA = cosm(g * a_dag_plus_a)
    
    c_dag_c_kin_L = sparse.kron(c_dag_c_kin_L, exp_iA)
    
    c_dag_c_kin_L += c_dag_c_kin_L.conjugate().T
    
    return c_dag_c_kin_L.tocsr()

def kin_L_red_rev (L, N, N_ph, g):
    
    states = states_gen(L , N) 
    #print(states)
    num_rows, num_cols = states.shape
    #print(states.shape)
    c_dag_c_kin_L = np.zeros((num_rows , num_rows))

    states_new = states.copy()
    #print(states_new)

    for i in range(num_rows): 
        #print('i is', i)
        for j in range(L - 1):
            #print('j is', j)
            if states_new[i][j] == 1 and states_new[i][j + 1] == 0:
                #print('old',states_new[i][:])
                states_new[i][j]     = 0
                states_new[i][j + 1] = 1
                #print('new', states_new[i][:])

                for k in range(num_rows):
                    if np.array_equal(states_new[i][:],states[k][:]):
                        #print('that state is', states_new[i][:])
                        #print('coordinates are', i , k)
                        c_dag_c_kin_L[i][k] = 1
                        break
                
                states_new[i][:] = states[i][:]
                #print('new to old', states_new[i][:])
                #print('old', states[i][:])

            else: 
                pass

    
    Mrx_a_dag   = Mrces_a(N_ph)[0]
    Mrx_a       = Mrces_a(N_ph)[1]
    
    a_dag_plus_a = Mrx_a_dag + Mrx_a
    
    exp_iA = expm(-1j * g * a_dag_plus_a)
    #exp_iA = cosm(g * a_dag_plus_a)
    
    c_dag_c_kin_L = sparse.kron(c_dag_c_kin_L, exp_iA)
    
    #c_dag_c_kin_L += c_dag_c_kin_L.conjugate().T
    
    return c_dag_c_kin_L.tocsr()

def WS_L_red (L, N, N_ph):
    
    states = states_gen(L , N) 
    #print(states)
    num_rows, num_cols = states.shape
    #print(states.shape)
    c_dag_c_WS_L = np.zeros((num_rows , num_rows))

    for i in range(num_rows): 
        #print('i is', i)
        cfnt = 0
        for j in range(L):
            #print('j is', j)
            if states[i][j] == 1 :
                cfnt  += (j + 1)
    
        c_dag_c_WS_L[i][i] = cfnt
        
    a_diag = np.eye(N_ph + 1)
    
    c_dag_c_WS_L = sparse.kron(c_dag_c_WS_L, a_diag)
            
    return c_dag_c_WS_L.tocsr()

def c_dag_c_i_red (L, N, N_ph, i):
    
    states = states_gen(L , N) 
    #print(states)
    num_rows, num_cols = states.shape
    #print(states.shape)
    c_dag_c_WS_L_i = np.zeros((num_rows , num_rows))

    for k in range(num_rows): 
        
        if states[k][i] == 1 :
            c_dag_c_WS_L_i[k][k] = 1
        
    a_diag = np.eye(N_ph + 1)
    
    c_dag_c_WS_L_i = sparse.kron(c_dag_c_WS_L_i, a_diag)
            
    return c_dag_c_WS_L_i.tocsr()




def Mrces_a(N_ph):
    a_dag_arr = np.zeros((N_ph + 1)**2)
    
    for i in range(N_ph):
        a_dag_arr[(i+1)* (N_ph + 1) + i] = np.sqrt(i + 1)

    Mrx_a_dag = np        . reshape   (a_dag_arr, (-1, N_ph + 1))
    Mrx_a     = Mrx_a_dag . transpose (                         )
        
    return Mrx_a_dag, Mrx_a



def c_dag_c_kin_L (g, N_ph, L, BC, OBC, PBC):
    
    c_dag_arr = np.array([[0, 0], [1, 0]])
    c_arr     = np.array([[0, 1], [0, 0]])
    
    I_JW      = np.array([[1, 0], [0, -1]])
    
    c_dag_c_kin = np.kron(c_dag_arr, c_arr)
    
    c_dag_c_kin_L = []
    
    for i in range(L): 
        
        
            
        Id_M_bf    = sparse.identity(2**i)
        Id_M_af    = sparse.identity(2**(L - i - 1))
        

        c_dag_c_i = sparse.kron(Id_M_bf, c_dag_c_kin)   
        c_dag_c_i = sparse.kron(c_dag_c_i, Id_M_af)
            
                
        c_dag_c_kin_L.append(c_dag_c_i) 
    
    if BC == PBC:
        
        I_JW_pbc = I_JW.copy()
        
        for j in range (L - 2):
            
            I_JW_pbc = sparse.kron(I_JW_pbc, I_JW) 
        
        c_dag_c_i = sparse.kron(c_dag_arr, I_JW_pbc)
        c_dag_c_i = sparse.kron(c_dag_c_i, c_arr   )
        
        c_dag_c_kin_L.append(c_dag_c_i) 
    
    else:
        pass
        
    Mrx   = (2**(L+1), 2**(L+1)) 
    A = sparse.csr_matrix(Mrx)
    #print(A.shape)   
    for i in c_dag_c_kin_L:
        #print(i.shape)   
            
        A = A + i
        
    Mrx_a_dag   = Mrces_a(N_ph)[0]
    Mrx_a       = Mrces_a(N_ph)[1]
    
    a_dag_plus_a = Mrx_a_dag + Mrx_a
    
    exp_iA = expm(1j * g * a_dag_plus_a)
    #exp_iA = cosm(g * a_dag_plus_a)
    
    A = sparse.kron(A, exp_iA)
    
    return A.tocsr()  


def c_c_dag_kin_L (g, N_ph, L, BC, OBC, PBC):
    
    c_dag_arr = np.array([[0, 0], [1, 0]])
    c_arr     = np.array([[0, 1], [0, 0]])
    
    I_JW      = np.array([[1, 0], [0, -1]])
    
    c_c_dag_kin = np.kron(c_arr, c_dag_arr)
    
    c_c_dag_kin_L = []
    
    for i in range(L): 
        
        
            
        Id_M_bf    = sparse.identity(2**i)
        Id_M_af    = sparse.identity(2**(L - i - 1))
        

        c_c_dag_i = sparse.kron(Id_M_bf, c_c_dag_kin)   
        c_c_dag_i = sparse.kron(c_c_dag_i, Id_M_af)
            
                
        c_c_dag_kin_L.append(c_c_dag_i) 
    
    if BC == PBC:
        
        I_JW_pbc = I_JW.copy()
        for j in range (L - 2):
            
            I_JW_pbc = sparse.kron(I_JW_pbc, I_JW) 
        
        c_c_dag_i = sparse.kron(c_arr, I_JW_pbc)
        c_c_dag_i = sparse.kron(c_c_dag_i, c_dag_arr   )
        
        c_c_dag_kin_L.append(c_c_dag_i)
    
    else:
        pass
        
    Mrx   = (2**(L+1), 2**(L+1)) 
    A = sparse.csr_matrix(Mrx)
        
    for i in c_c_dag_kin_L:
        
            
        A = A + i
        
    Mrx_a_dag   = Mrces_a(N_ph)[0]
    Mrx_a       = Mrces_a(N_ph)[1]
    
    a_dag_plus_a = Mrx_a_dag + Mrx_a
    
    exp_iA = expm(-1j * g * a_dag_plus_a)
    #exp_iA = cosm(g * a_dag_plus_a)

    
    A = sparse.kron(A, exp_iA)
    
    return A.tocsr()  


def c_dag_c_on_site_L (N_ph, L):
    
    c_dag_c_arr = np.array([[0, 0], [0, 1]])   
    #I_JW      = np.array([[1, 0], [0, -1]])
    

    
    c_dag_c_on_site_L = []
    
    for i in range(L + 1): 
            
        Id_M_bf    = sparse.identity(2**(i))
        Id_M_af    = sparse.identity(2**(L - i))
        

        c_dag_c_i = sparse.kron(Id_M_bf, c_dag_c_arr)   
        c_dag_c_i = sparse.kron(c_dag_c_i, Id_M_af) * (i + 1)

        c_dag_c_on_site_L .append(c_dag_c_i) 
        
    Mrx   = (2**(L + 1), 2** (L + 1)) 
    A = sparse.csr_matrix(Mrx)
    #print(A.shape)
        
    for i in c_dag_c_on_site_L:
        #print(i.shape)
        
        #print(i.shape)   
        A = A + i
        
    a_diag = np.eye(N_ph + 1)
    
    A = sparse.kron(A, a_diag)
        
    return A.tocsr()  


def a_dag_a (N_ph, L):
    
    Mrx_a_dag   = Mrces_a(N_ph)[0]
    Mrx_a       = Mrces_a(N_ph)[1]
    
    Mrx_a_dag_a = np.matmul(Mrx_a_dag, Mrx_a)
    
    Id_M_L   = sparse.identity(2**(L+1))
    
    a_dag_a = sparse.kron(Id_M_L, Mrx_a_dag_a)
    
    
    return a_dag_a.tocsr()

def a_dag_a_red (N_ph, L, N):
    
    Mrx_a_dag   = Mrces_a(N_ph)[0]
    Mrx_a       = Mrces_a(N_ph)[1]
    
    Mrx_a_dag_a = np.matmul(Mrx_a_dag, Mrx_a)
    
    dim      = np.int(math.factorial(L)/math.factorial(L-N)/math.factorial(N))
    Id_M_L   = sparse.identity(dim)
    
    a_dag_a = sparse.kron(Id_M_L, Mrx_a_dag_a)
    
    
    return a_dag_a.tocsr()

def a_dag_plus_a_red (N_ph, L, N):
    
    Mrx_a_dag   = Mrces_a(N_ph)[0]
    Mrx_a       = Mrces_a(N_ph)[1]
    
    a_dag_plus_a = Mrx_a_dag + Mrx_a
    dim =  np.int(math.factorial(L)/math.factorial(L-N)/math.factorial(N))
    Id_M_L   = sparse.identity(dim)
    
    a_dag_plus_a = sparse.kron(Id_M_L, a_dag_plus_a)
    
    return a_dag_plus_a.tocsr()


def c_dag_c_L (N_ph, L):
    
    c_dag_c_arr = np.array([[0, 0], [0, 1]])   
    #I_JW      = np.array([[1, 0], [0, -1]])
    

    
    c_dag_c_L = []
    
    for i in range(L + 1): 
        #print('i is', i)
            
        Id_M_bf    = sparse.identity(2**(i))
        Id_M_af    = sparse.identity(2**(L - i))
        

        c_dag_c_i = sparse.kron(Id_M_bf, c_dag_c_arr) 
        #print(c_dag_c_i)

        c_dag_c_i = sparse.kron(c_dag_c_i, Id_M_af)

        c_dag_c_L .append(c_dag_c_i) 
        
    Mrx   = (2**(L + 1), 2** (L + 1)) 
    A = sparse.csr_matrix(Mrx)
    #print(A.shape)
        
    for i in c_dag_c_L:
        #print(i.shape)
        
        #print(i.shape)   
        A = A + i
        
    a_diag = np.eye(N_ph + 1)
    
    A = sparse.kron(A, a_diag)
        
    return A.tocsr()  

def c_dag_c_i (N_ph, L, i):
    
    c_dag_c_arr = np.array([[0, 0], [0, 1]])   

            
    Id_M_bf    = sparse.identity(2**(i))
    Id_M_af    = sparse.identity(2**(L - i))
        

    c_dag_c_i = sparse.kron(Id_M_bf, c_dag_c_arr) 
    #print(c_dag_c_i)

    c_dag_c_i = sparse.kron(c_dag_c_i, Id_M_af)

    a_diag = np.eye(N_ph + 1)
    
    A = sparse.kron(c_dag_c_i, a_diag)
        
    return A.tocsr() 


def a_dag_plus_a (N_ph, L):
    
    Mrx_a_dag   = Mrces_a(N_ph)[0]
    Mrx_a       = Mrces_a(N_ph)[1]
    
    a_dag_plus_a = Mrx_a_dag + Mrx_a

    Id_M_L   = sparse.identity(2**(L+1))
    
    a_dag_plus_a = sparse.kron(Id_M_L, a_dag_plus_a)
    
    return a_dag_plus_a.tocsr()

def expectation_value (op, A):
    
    sA     = sparse.csr_matrix(A)
    oc2    = op.dot(sA.T)
    #oc2    = oc2.dot(sA.T)
    oc2    = np.conj(sA).dot(oc2)
    value= np.real_if_close(oc2[0,0])
    return value



class Vector:
    def __init__(self, A):
        self.v=sparse.csr_matrix(A)
    def expectation_value(self, op):
         oc2    = op.dot(self.v.T)
         oc2    = np.conj(self.v).dot(oc2)
         value= np.real_if_close(oc2[0,0])
         return value
    def apply(self, op):
        self.v=op.dot(self.v.T).transpose()

def time_evol (sA, U, del_t):
    
    #sA     = sparse.csr_matrix(A)
    #print(sA.shape)

    #U = sparse.linalg.expm(- H * del_t)
    #print(U.shape)
    sA_new = U.dot(sA.T)
    sA_new = sA_new.transpose()
    #A_new = sparse.tensordot(U,sA,axes=(1,0))
    return sA_new

# def a_avg_oc2 (A, N_ph, N_bn):
    
#     op     = a_oc_op(N_ph, N_bn)
#     sA     = sparse.csr_matrix(A)
#     oc2    = op.dot(op)
#     oc2    = oc2.dot(sA.T)
#     oc2    = np.conj(sA).dot(oc2)
    
#     return oc2[0,0]


#print(c_dag_c_kin_L(1,3,4))
#print(c_dag_c_on_site_L(3,10))

# BC = OBC

# h = 1.0
# g = 1.0
# Omega = 1.0 
# t_h = -1.0
# L = 7#index of last site (counting from 0)
# N_ph = 8

# k =  h * c_dag_c_L (N_ph, L)
# kk = c_c_dag_kin_L (g, N_ph, L, BC, OBC, PBC)
# kkk = c_dag_c_kin_L (g, N_ph, L, BC, OBC, PBC).tocsr() + k.tocsr()
# #kkk = c_c_dag_kin_L (g, N_ph, L, BC, OBC, PBC) + c_dag_c_kin_L (g, N_ph, L, BC, OBC, PBC) + h * c_dag_c_L (N_ph, L) 
# print(kkk)
# print('k is', k.getformat())
# print('kk is', kk.getformat())
# print('kkk is',kkk.getformat())
# g_arr = np.arange(0.0, 2.1, 0.1 )
# E0_lst = []
# a_oc_lst = []
# op = a_dag_a(N_ph, L)
# op1= c_dag_c_L(N_ph, L)
# eig_n = 40

# for g in g_arr:
#     print(g)
    
#     #H  = t_h * (c_dag_c_kin_L (g, N_ph, L, BC, OBC, PBC) + c_c_dag_kin_L (g, N_ph, L, BC, OBC, PBC)) + h * c_dag_c_on_site_L (N_ph, L) + Omega * a_dag_a(N_ph, L)
#     kin = t_h * c_dag_c_kin_L (g, N_ph, L, BC, OBC, PBC) 
#     kin_hc = kin.conjugate().T
#     H  = kin + kin_hc + h * c_dag_c_L (N_ph, L) + Omega * a_dag_a(N_ph, L)
#     #H  = kin + kin_hc + h * c_dag_c_on_site_L (N_ph, L) + Omega * a_dag_a(N_ph, L)

#     #print(H)

#     #eigenvalues, eigenvectors = eigsh(H, k = 10 , which = 'LA')
#     eigenvalues, eigenvectors = eigsh(H, k = eig_n , which = 'SA')
   
#     #E0_lst.append(eigenvalues[0])
    
#     # GS = eigenvectors[:,0]
#     # S9 = eigenvectors[:,9]
#     # a_oc = expectation_value ( op , GS )
#     # c_oc = expectation_value ( op1, GS)
#     for i in range (eig_n - 1):
#         S_i = eigenvectors[:,i]
#         c_oc = expectation_value ( op1, S_i)
#         print(c_oc)
#         eps = np.abs(c_oc - (L + 1)/2)
#         #print(eps)
#         if eps  < 10**(-4):
#             a_oc = expectation_value ( op , S_i )
#             a_oc_lst.append(a_oc)
#             print('filling is', c_oc)
#             break
#         else:
#             pass
        
#     #print(eigenvalues)    
    

# # plt.figure()

# # plt.plot(g_arr, E0_lst )
# # plt.xlabel(r'$g$')
# # plt.ylabel(r'$E_0$')
# # plt.legend()
# # #plt.savefig('E0_vs_t_pars1.png',  dpi = 1200)
# # plt.show()

# g_arr_GM = np.arange(0.0, 2.0, 0.2 )
# data_GM = np.load('N_avg_bos__TEBDPsi_GS_Nmax_8L_8g_1.8Omega_1J_1h_1V_0.npy')

# plt.figure()

# plt.plot(g_arr[1:], a_oc_lst[:])
# plt.plot(g_arr_GM, data_GM)
# plt.xlabel(r'$g$')
# plt.ylabel(r'$<N_{ph}>$')
# plt.legend()
# #plt.savefig('E0_vs_t_pars1.png',  dpi = 1200)
# plt.show()


