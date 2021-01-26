# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 08:59:27 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt

t=np.linspace(0, 2, 200)
x_t=np.load('C:/Users/giaco/Desktop/Cluster/Coherent_state/LC_coherent_L20_g0.25_Omega_10displacement_0.25dt_0.01X(t).npy')
plt.plot(t, x_t)
plt.xlabel('t')
plt.ylabel('<X(t)>')

"""
plt.text(0,0,r'$\Omega=10$')
plt.text(0,-0.10,r'$g=0.25$')
"""
plt.show()