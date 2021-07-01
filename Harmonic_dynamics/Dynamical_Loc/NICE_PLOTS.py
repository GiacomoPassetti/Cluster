# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 20:08:10 2021

@author: giaco
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
from scipy.linalg import expm, sinm, cosm, eigh
from scipy.sparse.linalg import eigsh
import matplotlib.patches as mpatches

from numpy import exp,arange
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from matplotlib import colors


import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm




def expectation(v, op):
    val = np.real_if_close(np.tensordot(v.conj().T, np.tensordot(op, v, 1),1))
    return val

def squeezednes_i(g, L, Omega, Nmax):
    B = np.diag(np.sqrt(np.arange((Nmax -1)) + 1), +1)
    Bd = np.diag(np.sqrt(np.arange((Nmax -1)) + 1), -1)
    Nb = Bd.dot(B) 
    X = B + Bd
    X= B+Bd
    Y= 1j*(Bd-B)
    XX=X.dot(X)
    YY=Y.dot(Y)    
    H=-cosm((g/np.sqrt(L))*X)*(L/np.pi)+(Omega*Nb)
    w, v= eigsh(H, 1, which='SA')
    v=v[:,0]
    dx=expectation(v, XX)
    dx=dx-(expectation(v, X)**2)

    dy=expectation(v, YY)

    dy=dy-(expectation(v, Y)**2)

    return np.sqrt(dy/dx)

def Squeezedness_plot(g_0,g_f,Omega1,Omega2,Omega3, L, Omega, Nmax):
    g=np.arange(g_0,g_f+0.025,.025 )
    Ls=[500]

    y1=[]
    for i in range(len(g)):
       y1.append(squeezednes_i(g[i], L, Omega1, Nmax))

 
    y01=[]
    for i in range(len(g)):
       y01.append(squeezednes_i(g[i], L, Omega2, Nmax))

 
    y10=[]
    for i in range(len(g)):
       y10.append(squeezednes_i(g[i], L, Omega3, Nmax))

    
    fig,ax = plt.subplots(dpi=600)
    gi=[0.2,0.75]
    ysq=[squeezednes_i(0.2, L, Omega1, Nmax), squeezednes_i(0.75, L, Omega1, Nmax)]
    fig,ax = plt.subplots(dpi=600)
    ax.plot(g, y01, color='lightskyblue', label=r'$\omega_{0}= 0.1 t_{h}$', zorder=1)
    ax.plot(g, y1, color='black', label=r'$\omega_{0}= t_{h}$', zorder=1)
    ax.plot(g, y10, color='peru', label=r'$\omega_{0}= 10 t_{h}$', zorder = 1)
    ax.scatter(gi,ysq,  s=50,marker='x',color='r',linewidths=1.5, zorder=2)
    return ax

def Ellypse_squeezing(g, L, Nmax, Omega):
    B = np.diag(np.sqrt(np.arange((Nmax -1)) + 1), +1)
    Bd = np.diag(np.sqrt(np.arange((Nmax -1)) + 1), -1)
    Nb = Bd.dot(B) 
    X = B + Bd
    X= B+Bd
    Y= 1j*(Bd-B)
    XX=X.dot(X)
    YY=Y.dot(Y)
    H=-cosm((g/np.sqrt(L))*X)*(L/np.pi)+(Omega*Nb)
    w, v= eigsh(H, 1, which='SA')
    v=v[:,0]
    dx=expectation(v, XX)
    dx=np.sqrt(dx-(expectation(v, X)**2))

    dy=expectation(v, YY)

    dy=np.sqrt(dy-(expectation(v, Y)**2))
    Ely = Ellipse((0, 0), dx, dy, color = 'black', fill = False)



    fig, ax = plt.subplots(dpi=800) 
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.arrow(-0.75*dx, 0, 1.5*dx, 0, head_width=0.05, head_length=0.05, fc='k', ec='k')
    ax.arrow(0, -0.75*dy, 0, 1.5*dy, head_width=0.05, head_length=0.05, fc='k', ec='k')
    ax.add_patch(Ely)
    return ax

def shaded_ellipse(g, L, Nmax, Omega):
    B = np.diag(np.sqrt(np.arange((Nmax -1)) + 1), +1)
    Bd = np.diag(np.sqrt(np.arange((Nmax -1)) + 1), -1)
    Nb = Bd.dot(B) 
    X = B + Bd
    X= B+Bd
    Y= 1j*(Bd-B)
    XX=X.dot(X)
    YY=Y.dot(Y)
    H=-cosm((g/np.sqrt(L))*X)*(L/np.pi)+(Omega*Nb)
    w, v= eigsh(H, 1, which='SA')
    v=v[:,0]
    dx=expectation(v, XX)
    dx=np.sqrt(dx-(expectation(v, X)**2))

    dy=expectation(v, YY)

    dy=np.sqrt(dy-(expectation(v, Y)**2))
    sigmax = dx
    sigmay = dy
    dx = 3*sigmax
    dy = 3*sigmay
    def z_func(x,y):

      z=(np.exp(-(x**2)/(2*sigmax**2)))*(np.exp(-(y**2)/(2*sigmay**2)))

      return z
    cnorm = colors.Normalize(-1, 1)
    x = arange(-dx,dx,dx*0.005)
    y = arange(-dy,dy,dy*0.005)
    X,Y = meshgrid(x, y) # grid of point
    Z = z_func(X, Y) # evaluation of the function on the grid
    
    fig, ax = plt.subplots(dpi = 600)
    ax.pcolormesh(X, Y, Z, norm=cnorm, cmap='RdBu')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    return ax
 



g = 1
L = 200
Omega1 = 0.1
Omega2 = 1
Omega3 = 10
g_0 = 0
g_f = 1
Nmax = 60
Omega = 1

fig = plt.subplots(dpi = 800)
#ax = Ellypse_squeezing(g, L, Nmax, Omega)
#ax = Squeezedness_plot(g_0,g_f,Omega1,Omega2,Omega3, L, Omega, Nmax)
#ax = shaded_ellipse(g, L, Nmax, Omega)



