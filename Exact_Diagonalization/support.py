# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 14:04:49 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
from scipy.linalg import expm, sinm, cosm, eigh
from tenpy.networks.site import FermionSite
from tenpy.networks.site import BosonSite
import copy


print(5)
t=np.arange(0,5, 0.01)
plt.plot(t,t,t,t**2,t,t**3)
plt.legend(['t', 'ts', 'tq'])
plt.show()