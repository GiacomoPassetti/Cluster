# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 11:03:06 2021

@author: giaco
"""

import matplotlib.pyplot as plt
import numpy as np

Omega=10
L=12
a=np.load('C:/Users/giaco/Desktop/Cluster/CDW/Quench_Exact/Quench/Data/Quench_Ui0.5_Uf_1.5dU_0.01imp_0.0Omega_10g_1.7dt_0.1E_t.npy')
t=np.arange(0,10,0.001)
plt.plot(t,a)
plt.xlabe('t')
plt.ylabel('E(t)-E_gs(U(t))')


