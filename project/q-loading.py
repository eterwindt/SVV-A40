# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:05:27 2020

@author: sebas
"""

#----------------------------- imports ----------------------------------------
from numpy import pi, cos

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
#----------------------- functions definitions --------------------------------
Ca = 0.505 #[m], chord of aileron
la = 1.611 #[m], span of aileron

#-------------------------- main program --------------------------------------
filename = 'aerodynamicloadf100.dat'

def readfile(filename):
    file = open(filename, "r")
    lines = file.readlines()                    # read all lines, returns list
    for n in range(0,len(lines)):
        lines[n] = lines[n].split(",")          # Split lines at "," to get individual entries
        lines[n][-1] = lines[n][-1][:-1]        # Remove \n at last entries
    file.close()
    return lines

q = readfile(filename)                          # read aerodynamic load, list of 81 rows, 41 collumns
Nx = len(q[0])                                  # Amount of entries in x-direction, 41 collumns
Nz = len(q)                                     # Amount of entries in z-direction, 81 rows

theta_z = []
z       = []
theta_x = []
x       = []
for i in range(1,Nz+1):                         # Get a list of the z coordinate of every entry
    theta_z_i = (i-1)*pi/Nz
    theta_z_i1 = (i)*pi/Nz

    z_i = - 0.5*((Ca/2)*(1-cos(theta_z_i)) + (Ca/2)*(1-cos(theta_z_i1)))
    z.append(z_i)

for j in range(1,Nx+1):                         # Get a list of the x coordinate of every entry
    theta_x_j = (j-1)*pi/Nx
    theta_x_j1 = (j)*pi/Nx

    x_j = 0.5*((la/2)*(1-cos(theta_x_j)) + (la/2)*(1-cos(theta_x_j1)))
    x.append(x_j)

#------------------------ 3d surface plot -------------------------------------

X = np.float64(x)
Y = np.float64(z)
Z = np.float64(q)
X2 = np.vstack([X]*Nz)
Y2 = np.transpose(np.vstack([Y]*Nx))

'''Mesh grid to see where the Aerodynamic Loads point are on the aileron'''
fig = plt.figure()
x1, y1 = np.meshgrid(X,Y)
plt.plot(x1,y1, marker='.', markersize=2, color='k', linestyle='none')
plt.show()

''' 3D Surface plot of aerodynamic load with 1:1 axis '''
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X2, Y2, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set_xlim(0, 1.7)
ax.set_ylim(-1, 0.7)
ax.set_xlabel('X')
ax.set_ylabel('Z')
ax.set_zlabel('q')
fig.colorbar(surf, shrink=0.5, aspect=5)



