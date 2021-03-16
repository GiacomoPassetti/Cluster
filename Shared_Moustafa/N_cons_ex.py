# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 10:31:01 2021

@author: giaco
"""

import numpy as np
import itertools
import math
from scipy import sparse
from scipy.sparse.linalg import expm
from scipy.sparse.linalg import eigsh 
import matplotlib.pyplot as plt 
from scipy.linalg import cosm

def states_gen(L,N):
    which = np.array(list(itertools.combinations(range(L), N)))
    #print(which)
    grid = np.zeros((len(which), L), dtype="int8")

    # Magic
    grid[np.arange(len(which))[None].T, which] = 1
    
    return grid


def kin_L_red (L, N, N_ph, g):
    
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

def dN_dN (L, N, N_ph, g):
    
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
            if states_new[i][j] == 1 and states_new[i][j + 1] == 1:
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

def Mrces_a(N_ph):
    a_dag_arr = np.zeros((N_ph + 1)**2)
    
    for i in range(N_ph):
        a_dag_arr[(i+1)* (N_ph + 1) + i] = np.sqrt(i + 1)

    Mrx_a_dag = np        . reshape   (a_dag_arr, (-1, N_ph + 1))
    Mrx_a     = Mrx_a_dag . transpose (                         )
        
    return Mrx_a_dag, Mrx_a