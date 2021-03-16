# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 14:23:26 2021

@author: giaco
"""

import tenpy
import copy
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/CDW')
from N_cons import Suz_trot_im, H_Peier_bond, product_state, Iterative_g, Iterative_g_from_load
import numpy as np
import numpy.linalg as alg
import matplotlib.pyplot as plt
from tenpy import models
from tenpy.networks.site import SpinSite
from tenpy.networks.site import FermionSite
from tenpy.networks.site import BosonSite
from tenpy.models.model import CouplingModel
from tenpy.models.model import CouplingMPOModel
from tenpy.models.xxz_chain import XXZChain
from tenpy.models.spins import SpinModel
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS
from tenpy.models.lattice import Lattice
from tenpy.tools.params import get_parameter
import tenpy.linalg.np_conserved as npc
from tenpy.networks.mpo import MPO, MPOEnvironment
import tenpy.linalg.charges as charges
from tenpy.models.lattice import Chain
from scipy.linalg import expm
from tenpy.models.fermions_spinless import FermionModel
from tenpy.algorithms.tebd import Engine
import pickle
from tenpy.linalg.charges import LegCharge, ChargeInfo
from tenpy.algorithms.tebd import Engine

def example_DMRG_heisenberg_xxz_infinite(Jz, L,  ):
    print("finite DMRG, Heisenberg XXZ chain")
    
    model_params = dict(
        L=L,
        S=0.5,  # spin 1/2
        Jxx=1.,
        
        Jz=Jz,  # couplings
        bc_MPS='finite',
        )
    M = XXZChain(model_params)
    
    product_state = ["up", "down"]*(L//2)  # initial Neel state
    psi = MPS.from_product_state(M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)
    dmrg_params = {
        'mixer': True,  # setting this to True helps to escape local minima
        'trunc_params': {
            'chi_max': 100,
            'svd_min': 1.e-10,
        },
        'max_E_err': 1.e-10,
    }
    info = dmrg.run(psi, M, dmrg_params)
    E = info['E']
    print("E = {E:.13f}".format(E=E))
    print("final bond dimensions: ", psi.chi)
    Sz = psi.expectation_value("Sz")  # Sz instead of Sigma z: spin-1/2 operators!
    mag_z = np.mean(Sz)
    print("<S_z> = [{Sz0:.5f}, {Sz1:.5f}]; mean ={mag_z:.5f}".format(Sz0=Sz[0],
                                                                     Sz1=Sz[1],
                                                                     mag_z=mag_z))
    # note: it's clear that mean(<Sz>) is 0: the model has Sz conservation!
    
    corrs = psi.correlation_function("Sz", "Sz", sites1=[0],sites2=range(10))
    print("correlations <Sz_i Sz_j> =")
    print(corrs)
    return E, psi, M, corrs


E, psi, M, corrs = example_DMRG_heisenberg_xxz_infinite(Jz=0, L=20)
plt.plot(psi.expectation_value('Sz'))
