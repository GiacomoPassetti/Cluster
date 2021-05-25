# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 10:17:49 2021

@author: giaco
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quantum_information/Spin_v')

epsilon=0.01
J=1.0
tmax=10
steps=5
dt=0.01
ts=np.arange(0,10+steps*dt, steps*dt)
hmax=20.0
L=20

ts=np.arange(0,tmax, 2*dt*steps)
#a=np.load("C:/Users/giaco/Desktop/Cluster/Quantum_information/Data/Random_QI_Algorithm_dt0.01t_max10epsilon_"+str(epsilon)+"h_max"+str(hmax)+"chi_max300steps_5L_"+str(L)+"SMax.npy")
#b=np.load("C:/Users/giaco/Desktop/Cluster/Quantum_information/Data/Random_QI_Algorithm_dt0.01t_max10epsilon_"+str(epsilon)+"h_max"+str(hmax)+"chi_max300steps_5L_"+str(L)+"Smean.npy")
a=np.load("C:/Users/giaco/Desktop/Cluster/Quantum_information/Spin_v/Random_QI_Algorithm_dt0.01t_max10epsilon_0.01h_max1chi_max300steps_10L_20err.npy")
b=np.load("C:/Users/giaco/Desktop/Cluster/Quantum_information/Spin_v/Random_QI_Algorithm_dt0.01t_max10epsilon_0.01h_max5chi_max300steps_10L_20err.npy")
c=np.load("C:/Users/giaco/Desktop/Cluster/Quantum_information/Spin_v/Random_QI_Algorithm_dt0.01t_max10epsilon_0.01h_max10chi_max300steps_10L_20err.npy")
plt.plot(ts, a)
#plt.plot(ts, b)
#plt.plot(ts, c)




