# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 18:58:32 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt

t=np.linspace(0, 2, 200)
x_t=np.load('LC_coherent_L20_g0_Omega_10displacement_1dt_0.005X(t)')
plt.plot(t, x_t)
plt.xlabel('t')
plt.ylabel('<X(t)>')

"""
plt.text(0,0,r'$\Omega=10$')
plt.text(0,-0.10,r'$g=0.25$')
"""
plt.show()