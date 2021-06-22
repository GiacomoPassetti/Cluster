# -*- coding: utf-8 -*-
"""
Created on Tue May 25 08:16:02 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt
import sys


ts = np.arange(0,1.025, 0.025)
Fold = 'C:/Users/giaco/Desktop/Cluster/Quantum_information/Entropy_propagation/Data_ent_lev/'


a = "XXZ_spin_propagation_Delta_0.5.npy"
    


data= np.load(Fold+a)


print(data)
    