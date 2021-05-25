# -*- coding: utf-8 -*-
"""
Created on Mon May  3 15:49:54 2021

@author: giaco
"""
import numpy as np
import matplotlib.pyplot as plt 


fig,ax = plt.subplots(dpi=600)


# Have a look at the colormaps here and decide which one you'd like:
# http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html








fig,ax = plt.subplots(dpi=600)
x= np.linspace(-np.pi, np.pi, 200)
y1 =- np.cos(x)
#y1 = np.sin(x)

line1= plt.axvline((-np.pi/2),-1, 0.95, color="red")
line2= plt.axvline((np.pi/2), -1, 0.95, color="red")
line3=plt.hlines(0.99, (-np.pi/2), (np.pi/2), color="red")

ax.plot(x, y1, color='black')
#ax.plot(x, y2, color='black')







    
    








