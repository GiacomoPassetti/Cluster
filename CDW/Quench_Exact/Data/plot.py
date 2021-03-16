# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 09:01:37 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/CDW/Quench_Exact/Data')

a=np.load('C:/Users/giaco/Desktop/Cluster/CDW/Quench_Exact/Quench/E_t.npy')

plt.plot(a)


