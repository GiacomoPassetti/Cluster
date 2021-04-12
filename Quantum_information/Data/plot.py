# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 10:17:49 2021

@author: giaco
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quantum_information/Data')

epsilon=0.1
J=1.0
tmax=10
steps=5
dt=0.01
ts=np.arange(0,10+steps*dt, steps*dt)
hmax=20.0
L=20

"""
a=np.load("C:/Users/giaco/Desktop/Cluster/Quantum_information/Data/Mean_entropy_dt0.01t_max10epsilon_"+str(epsilon)+"h_max"+str(hmax)+"chi_max300steps_5J_"+str(J)+"SMax.npy")
b=np.load("C:/Users/giaco/Desktop/Cluster/Quantum_information/Data/Mean_entropy_dt0.01t_max10epsilon_"+str(epsilon)+"h_max"+str(hmax)+"chi_max300steps_5J_"+str(J)+"Smean.npy")
"""


"First plot: J=0 and increasing the the error in the swap gate:"

"""
f=[]
for epsilon in [0.0 , 0.0001, 0.001, 0.1]:
    f.append(np.load("C:/Users/giaco/Desktop/Cluster/Quantum_information/Data/Mean_entropy_dt0.01t_max10epsilon_"+str(epsilon)+"h_max"+str(hmax)+"chi_max300steps_5J_"+str(J)+"Smean.npy"))
    
plt.plot(ts,f[0], ts, f[1], ts, f[2], ts, f[3])
plt.xlabel("t")
plt.ylabel(r"$S_{avg}$")
plt.legend([r"$\epsilon=0$", r"$\epsilon=0.0001$",r"$\epsilon=0.001$", r"$\epsilon=0.1$"])
plt.title(r"$|h_{max}|=5$")
"""

"""
aa=[]
for epsilon in [0.1, 0.01, 0.001, 0.0001]:
 avgs=[]
 for RM in range(10):
    avgs.append(np.load("C:/Users/giaco/Desktop/Cluster/Quantum_information/Data/Mean_entropy_dt0.01t_max10epsilon_"+str(epsilon)+"h_max"+str(hmax)+"chi_max300steps_5J_"+str(J)+"Random_run"+str(1+RM)+"Smean.npy"))
 avged=[]

 for i in range(201):
    a0=0
    for j in range(10):
        a0=a0+avgs[j][i]
        a0=a0/10
        
    avged.append(a0)
 aa.append(avged)
plt.figure(dpi=1200)
plt.plot(ts, aa[0], ts, aa[1], ts, aa[2], ts, aa[3])
plt.xlabel(r"$t$")
plt.ylabel(r"$S_{mean}$")
plt.title("Average Entropy 40 sites  |h|=10")
plt.legend([r"$\epsilon=0.1$", r"$\epsilon=0.01$", r"$\epsilon=0.001$", r"$\epsilon=0.0001$"])
"""

aa=[]
for hmax in [10.0, 20.0]:
 avgs=[]
 for RM in range(10):
    avgs.append(np.load("C:/Users/giaco/Desktop/Cluster/Quantum_information/Data/Mean_entropy_dt0.01t_max10epsilon_"+str(epsilon)+"h_max"+str(hmax)+"chi_max300steps_5J_"+str(J)+"Random_run"+str(1+RM)+"Smean.npy"))
 avged=[]

 for i in range(201):
    a0=0
    for j in range(10):
        a0=a0+avgs[j][i]
        a0=a0/10
        
    avged.append(a0)
 aa.append(avged)
plt.figure(dpi=1200)
plt.plot(ts, aa[0], ts, aa[1])
plt.xlabel(r"$t$")
plt.ylabel(r"$S_{mean}$")
plt.title("Average Entropy 40 sites  eps=0.1")
plt.legend([r"$|h|=5$", r"|h|=10"])