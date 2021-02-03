# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 18:58:32 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt

X=[np.load()]
t=np.arange(0,0.5,0.1)
plt.plot(t,X[0],t,X[1],t,X[2])
plt.xlabel('t')
plt.ylabel(r'$<X(t)>    \Omega=$'+str(Omega))
plt.legend(['g=0.25', 'g=0.75', 'g=1.5'])