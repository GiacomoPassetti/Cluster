# -*- coding: utf-8 -*-
"""
Created on Thu May 27 15:18:08 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt

R=0.8

def I(x, R):
    X= np.sqrt(R)
    I = 1/(((1-X)**2) + 4*X*(np.sin(x))**2)
    return I

x = np.linspace(2,10.5, 800)
y = I(x, R)


fig, ax = plt.subplots(dpi=1000) 
# note we must use plt.subplots, not plt.subplot
# (or if you have an existing figure)
# fig = plt.gcf()
# ax = fig.gca()

ax.vlines(np.pi, 0, 90,"cornflowerblue", linewidth = 3)
ax.text(np.pi+0.1, 80, r'$\omega_{0}$', fontsize=20)
ax.plot(x,y, "black", linewidth = 0.9)
ax.set_aspect(0.05)
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)





