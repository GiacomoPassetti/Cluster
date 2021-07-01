# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 15:50:30 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('C:/Users/giaco/Desktop/Cluster/Quantum_information/Err_sx/Data')

eps = 0.01
iterations = 4
avg1 = []
avg2 = []
avg3 = []
for i in range(iterations):                                                              #Random_algorithm_J_Error_sx_1eps_0.01Iterations_4t_max3Random_IT_1hMAX_0.2avg_1
      avg1.append(np.load('C:/Users/giaco/Desktop/Cluster/Quantum_information/Err_sx/Data/Random_algorithm_J_Error_sx_1eps_'+str(eps)+'Iterations_4t_max3Random_IT_'+str(i+1)+'hMAX_0.2avg_1.npy'))
      avg2.append(np.load('C:/Users/giaco/Desktop/Cluster/Quantum_information/Err_sx/Data/Random_algorithm_J_Error_sx_1eps_'+str(eps)+'Iterations_4t_max3Random_IT_'+str(i+1)+'hMAX_2.0avg_1.npy'))
      avg3.append(np.load('C:/Users/giaco/Desktop/Cluster/Quantum_information/Err_sx/Data/Random_algorithm_J_Error_sx_1eps_'+str(eps)+'Iterations_4t_max3Random_IT_'+str(i+1)+'hMAX_20.0avg_1.npy'))

err1 = zip(avg1[0], avg1[1])
err2 = zip(avg2[0], avg2[1])
err3 = zip(avg3[0], avg3[1])
err1 = [x + y for (x, y) in err1] 
err2 = [x + y for (x, y) in err2]
err3 = [x + y for (x, y) in err3]
for i in range(iterations - 2):
    err1 = zip(err1, avg1[2+i])
    err2 = zip(err2, avg2[2+i])
    err3 = zip(err3, avg3[2+i])
    err1 = [x + y for (x, y) in err1] 
    err2 = [x + y for (x, y) in err2]
    err3 = [x + y for (x, y) in err3]
    
 
err1 = [x/iterations for x in err1]
err2 = [x/iterations for x in err2]
err3 = [x/iterations for x in err3]


plt.figure(dpi=900)
t = np.arange(0,3.025,0.025)
plt.plot(t, err1, t, err2, t, err3)
plt.xlabel(r"$t$")
plt.ylabel(r"$|<\psi^{ex}(t)|\psi^{err}(t)>|$")
plt.legend([r"$hmax = 0.1$", r"$hmax = 1$", r"$hmax = 10$"])
plt.title(r"$\epsilon = $"+str(eps))


