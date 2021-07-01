# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 11:20:59 2021

@author: giaco
"""

import matplotlib.pyplot as plt
import numpy as np
from tenpy.networks.site import BosonSite
from scipy.linalg import expm, sinm, cosm, eigh
from scipy.sparse.linalg import eigsh
import matplotlib.patches as patches
import matplotlib.ticker as ticker

"""
#Square 1
delta = 0.04
deltav = 0.01
ticks = (-np.pi, -np.pi/2, 0, np.pi/2, np.pi)
fig, ax = plt.subplots(dpi = 800)

ax.set_xticks(ticks, minor=False)
ax.set_yticks((0,1), minor=False)

ax.set_xticklabels((r"$-\pi$", r"$-\frac{\pi}{2}$", r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"), fontsize=8)
ax.set_yticklabels((r"$0$", r"$1$"), fontsize=8)
ax.set_xlabel(r"$k$", fontsize = 8)
ax.set_ylabel(r"$n(k)$", fontsize=8)
fig.set_size_inches(2., 1.5)

plt.vlines(0, 0, 1.1, color = "grey", linestyle="--", linewidth=0.8)
arrow = patches.FancyArrowPatch((0, 1.1), (0.4, 1.1), arrowstyle='->', mutation_scale=10, zorder = 100, linewidth=1.5, color = 'black')
ax.add_patch(arrow)
plt.vlines(-np.pi/2, 0, 1+deltav, color = "black")
plt.vlines(np.pi/2, 0, 1+deltav, color = "black")
plt.hlines(0, -np.pi, -np.pi/2+delta, color = "black")
plt.hlines(1, -np.pi/2, np.pi/2, color = "black")
plt.hlines(0, np.pi/2-delta, np.pi, color = "black")
plt.ylim([-0.05, 1.5])
plt.text(-1.5, 1.3, "FS shift = 0", fontsize=8)

plt.show()
"""

"""
#Square 2
delta = 0.04
deltav = 0.01
ticks = (-np.pi, -np.pi/2, 0, np.pi/2, np.pi)
fig, ax = plt.subplots(dpi = 800)

ax.set_xticks(ticks, minor=False)
ax.set_yticks((0,1), minor=False)

ax.set_xticklabels((r"$-\pi$", r"$-\frac{\pi}{2}$", r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"), fontsize=8)
ax.set_yticklabels((r"$0$", r"$1$"), fontsize=8)
ax.set_xlabel(r"$k$", fontsize = 8)
ax.set_ylabel(r"$n(k)$", fontsize=8)
fig.set_size_inches(2., 1.5)

plt.vlines(np.pi/2, 0, 1.1, color = "grey", linestyle="--", linewidth=0.8)
arrow = patches.FancyArrowPatch((0, 1.1), (np.pi/2, 1.1), arrowstyle='->', mutation_scale=10, zorder = 100, linewidth=1.5, color = 'black')
ax.add_patch(arrow)
plt.vlines(0, 0, 1+deltav, color = "black")
plt.vlines(np.pi, 0, 1+deltav, color = "black")
plt.hlines(0, -np.pi, delta, color = "black")
plt.hlines(1, 0, np.pi, color = "black")

plt.ylim([-0.05, 1.5])
plt.text(-1.5, 1.3, r"FS shift = $\frac{\pi}{2}$", fontsize=8)

plt.show()
"""


#Square 3
delta = 0.04
deltav = 0.01
ticks = (-np.pi, -np.pi/2, 0, np.pi/2, np.pi)
fig, ax = plt.subplots(dpi = 800)

ax.set_xticks(ticks, minor=False)
ax.set_yticks((0,1), minor=False)

ax.set_xticklabels((r"$-\pi$", r"$-\frac{\pi}{2}$", r"$0$", r"$\frac{\pi}{2}$", r"$\pi$"), fontsize=8)
ax.set_yticklabels((r"$0$", r"$1$"), fontsize=8)
ax.set_xlabel(r"$k$", fontsize = 8)
ax.set_ylabel(r"$n(k)$", fontsize=8)
fig.set_size_inches(2., 1.5)

plt.vlines(np.pi, 0, 1.1, color = "grey", linestyle="--", linewidth=0.8)
arrow = patches.FancyArrowPatch((0, 1.1), (np.pi, 1.1), arrowstyle='->', mutation_scale=10, zorder = 100, linewidth=1.5, color = 'black')
ax.add_patch(arrow)
plt.vlines(-np.pi/2, 0, 1+deltav, color = "black")
plt.vlines(np.pi/2, 0, 1+deltav, color = "black")
plt.hlines(1, -np.pi, -np.pi/2, color = "black")
plt.hlines(1, np.pi/2, np.pi, color = "black")
plt.hlines(0, -np.pi/2-delta, np.pi/2+delta, color = "black")

plt.ylim([-0.05, 1.5])
plt.text(-1.5, 1.3, r"FS shift = $\pi$", fontsize=8)

plt.show()
