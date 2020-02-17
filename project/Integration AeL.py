# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:56:53 2020

@author: esmee
"""
#----------------------------- imports ----------------------------------------
from numpy import pi, cos
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

q = readfile(filename)                          #Importing all aerodynamic loading data into lists in lists
Nx = len(q[0])
Nz = len(q)

theta_z = []
z       = []
theta_x = []
x       = []
for i in range(1,Nz+1):
    theta_z_i = (i-1)*pi/Nz
    theta_z_i1 = (i)*pi/Nz

    z_i = - 0.5*((Ca/2)*(1-cos(theta_z_i)) + (Ca/2)*(1-cos(theta_z_i1)))
    z.append(z_i)

for j in range(1,Nx+1):
    theta_x_j = (j-1)*pi/Nx
    theta_x_j1 = (j)*pi/Nx

    x_j = 0.5*((la/2)*(1-cos(theta_x_j)) + (la/2)*(1-cos(theta_x_j1)))
    x.append(x_j)

data = np.array([np.array(xi) for xi in q])     #Coverting q into array

""" Now to first interpolate wrt to x"""

L0x = []

for p in range(len(x)):
    i = 1.
    while i < 40.:
        L0xp = np.prod(((float(x[p]) - float(data[0,i]))/(float(data[0,0]) - float(data[0,i]))))
        i = i+1.
    L0x.append(L0xp)

print(i)



