# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 11:19:41 2021

@author: giaco
"""

import numpy as np
import matplotlib.pyplot as plt

L = 30
hmax = 0.0
IT = 9
ts = np.arange(0,20, 0.025)
a = np.zeros(len(ts))

def List_merger(avg, iterations):
    err1 = zip(avg[0], avg[1])
    err1 = [x + y for (x, y) in err1] 
    for i in range(iterations - 2):
      err1 = zip(err1, avg[2+i])

    err1 = [x + y for (x, y) in err1] 
    err1 = [x/iterations for x in err1]
    return err1

    
 

    

dats = []
for i in[1,2,3,4,5,6,7,9,10]:
    dats.append(np.load('C:/Users/giaco/Desktop/MBL/MBL_transition/Data/MBLtrs_hmax'+str(hmax)+'L_30RI_'+str(i)+'.npy'))
ABent=[]
for i in range(len(dats)):
    for j in range(len(dats[i])):
        ABent.append(dats[i][j][int(L/2)])

        
