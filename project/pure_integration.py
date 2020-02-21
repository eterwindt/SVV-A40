# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 17:40:58 2020

@author: esmee
"""
import numpy as np


def g(x2):
    integr = []
    for j in range(len(x)-1):
        dothing = True
        integrsmall = []
        for i in range(len(z)-1):
            G = coefarray[j,i,0]*(z[i+1]-z[i])+coefarray[i,i,1]*x2*(z[i+1]-z[i])\
                +0.5*coefarray[j,i,2]*(z[i+1]**2-z[i]**2)+0.5*coefarray[i,i,3]*x2*(z[i+1]**2-z[i]**2)
            
            integrsmall.append(G)
            if dothing == True:
                dothing = False
        integrsmall = np.array(integrsmall)
#        integrsmall = np.reshape(integrsmall, (40,1))
        integrsmall = np.sum(integrsmall)
        
        integr.append(integrsmall)
    integr = np.array(integr)
    #integr = np.reshape(integr, (80,40))
    return integr


if __name__ == '__main__':
    x_test = 1.4
    test_int = g(x_test)