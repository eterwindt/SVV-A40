# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:56:53 2020

@author: esmee
"""
#----------------------------- imports ----------------------------------------
from numpy import pi, cos
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
from scipy.integrate import simps, trapz
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

q = np.array(q, dtype='f8', copy=True, order='K')


def coefarrayFunc(q, Ca, la):
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

    data = np.array(q, dtype='f8', copy=True, order='K')     #Coverting q into array
    
    z = np.array(z)*-1
    #x = np.linspace(0, 40, 41)
    #z = np.linspace(0, 80, 81)
    
    """ Now interpolate wrt to x and z"""

    coefarray = []
    for j in range(len(z)-1):
        dothing = True
        coefsmall = []
        for i in range(len(x)-1):
            
            intarray = [[1, float(x[i]),   float(z[j]),   float(x[i])*float(z[j])],         #bivariate linear spline interpolation matrix
                        [1, float(x[i]),   float(z[j+1]), float(x[i])*float(z[j+1])],
                        [1, float(x[i+1]), float(z[j]),   float(x[i+1])*float(z[j])],
                        [1, float(x[i+1]), float(z[j+1]), float(x[i+1])*float(z[j+1])]]

            dataarray = [[float(data[j,i])],                    #solution matrix from aerodynnamic loading
                        [float(data[j+1,i])],
                        [float(data[j,i+1])],
                        [float(data[j+1,i+1])]]
            
            coef = np.dot(np.linalg.inv(intarray),dataarray)    #obtaining coefficient matrix for solving linear splines
            
            coefsmall.append(coef)                  #creating a matrix with coefficients for all the data points on the aileron in x-dir
            if dothing == True:
                #print(coef)
                dothing = False
          
        
        coefsmall = np.array(coefsmall)
        coefsmall = np.reshape(coefsmall, (40, 1, 4, 1))
        
        coefarray.append(coefsmall)
    coefarray = np.array(coefarray, dtype='f8')                 #creating the matrix in z-dir by appending the x-dir matrices
    coefarray = np.reshape(coefarray, (80, 40,  4))

    return coefarray, x, z 

coefarray, x, z = coefarrayFunc(q, Ca, la)

"""
Numerical Verification of the linear bivariate spline interpolation by 
comparing the outcome of the function with the original data           
for j in range(len(z)-1):
        for i in range(len(x)-1):
            coefs = coefarray[j,i,:]           
            function = coefs[0] + coefs[1]*x[i] + coefs[2]*z[j] + coefs[3]*x[i]*z[j]
#            print (function)
"""     
        
"""Now to integrate the function wrt y to obtain a 1D f(x)"""
'''
def print_menu():
    print( "Choose an option")
    print() 
    print( "1) Find Rforce at an x-coordinate due to aerodynamic loading")
    print( "2) Force f(x) integrated once")
    print( "3) Force f(x) integrated twice")
    print( "4) Force f(x) integrated thrice")
    print( "5) Force f(x) integrated four times") '''
    

'''
def g(x1):
    integr = []
    
    for j in range(len(x)-1):
        dothing = True
        integrsmall = []
        for i in range(len(z)-1):
            a0 = float(coefarray[i,j,0])
            a1 = float(coefarray[i,j,1])
            a2 = float(0.5*coefarray[i,j,2])
            a3 = float(0.5*coefarray[i,j,3])
            
            part1 = a0*(z[i+1]-z[i])
            part2 = a1*(z[i+1]-z[i])
            part3 = a2*(z[i+1]**2-z[i]**2)
            part4 = a3*(z[i+1]**2-z[i]**2)
            
            if choice == 1:
                G = part1 + (x1)*part2 + part3 + (x1)*part4
                integrsmall.append(G)
                #Integrated function of aerodynamic loading, wrt z                
            elif choice == 2:
                G = part1*x1 + 0.5*part2*(x1**2) + part3*x1 + 0.5*part4*(x1**2)
                integrsmall.append(G)
                #Integrated function of x wrt x
            elif choice == 3:
                G = 0.5*part1*(x1**2) + (1/6)*part2*(x1**3) + 0.5*part3*(x1**2) + (1/6)*part4*(x1**3)
                integrsmall.append(G)
                #Integrated function of x wrt x twice
            elif choice == 4:
                G = (1/6)*part1*(x1**3) + (1/24)*part2*(x1**4) + (1/6)*part3*(x1**3) + (1/24)*part4*(x1**4)
                integrsmall.append(G)
                #Integrated function of x wrt x thrice
            elif choice == 5:
                G = (1/24)*part1*(x1**4) + (1/120)*part2*(x1**5) + (1/24)*part3*(x1**4) + (1/120)*part4*(x1**5)
                integrsmall.append(G)
                #Integrated function of x wrt x four times
            else: 
                print("That is not a valid choice")
                print("Too bad, try again")
                break
            if dothing == True:
                dothing = False
        integrsmall = np.array(integrsmall)
        integrsmall = np.sum(integrsmall) #summing the parts in the z direction to gain the force for a coordinate x
        
        integr.append(integrsmall)
       
    integr = np.array(integr, dtype = 'f8') #array containing integrated functions to get g(x)
    
    for i in range(len(x)-1):
        if float(value) > float(x[i]):       #determining in which part of the aileron we are and which g(x) to use
            pass
        else:
            t = g(value)
            finalg = t[i]
            print(finalg)  
    return integr
'''


def integrFunc(coord, coefarray, x, z, choice):
    #INTEGRATE WRT Z
    
    coeffarray2 = np.zeros((80, 40, 5))
    sumert1 = 0
    sumert2 = 0 
    xlist1 = []
    xlist2 = []
    xlist3 = []
    
    for j in range((max(np.argwhere(np.array(x) < coord)))[0]):
        zlist = [] 
        for i in range(len(z)-1):
            
            
            a0 = float(coefarray[i,j,0])
            a1 = float(coefarray[i,j,1])
            a2 = float(coefarray[i,j,2])
            a3 = float(coefarray[i,j,3])

            
            part1 = a0*(z[i+1]-z[i])
            part2 = a1*(z[i+1]-z[i])
            part3 = a2*(z[i+1]**2-z[i]**2)*1/2
            part4 = a3*(z[i+1]**2-z[i]**2)*1/2
            

            start = a0*z[i] + a1*z[i]*x[j] + a2*(z[i]**2)*(1/2) + a3*(z[i]**2)*(1/2)*x[j]
            end = a0*z[i+1] + a1*z[i+1]*x[j] + a2*(z[i+1]**2)*(1/2) + a3*(z[i+1]**2)*(1/2)*x[j] 

            if i == 0:
                const = 0 
            if i != 0:
                const = zlist[-1] - start
            start += const
            end += const
            
            zlist.append(end)
            
       

        
            
        #if j == 0:
            #G = part1 + (x1)*part2 + part3 + (x1)*part4 + const
        
            #integrsmall.append(G1)
            #Integrated function of aerodynamic loading, wrt z                
        
        startx1 = a0*z[-1]*x[j] + a1*0.5*z[-1]*(x[j]**2) + (0.5*a2*z[-1]**2)*x[j] + 0.25*a3*(z[-1]**2)*(x[j]**2) + const*x[j]
        endx1 = a0*z[-1]*x[j+1] + a1*0.5*z[-1]*(x[j+1]**2) + (0.5*a2*z[-1]**2)*x[j+1] + 0.25*a3*(z[-1]**2)*(x[j+1]**2) + const*x[j+1]
        
        if j == 0:
            constx1 = 0 
        if j != 0:
            constx1 = xlist1[-1] - startx1
        startx1 += constx1 
        endx1 += constx1

        xlist1.append(endx1)
        
        startx2 = (1/2)*a0*z[-1]*(x[j]**2) + (1/3)*a1*0.5*z[-1]*(x[j]**3) + (0.25*a2*z[-1]**2)*(x[j]**2) + (1/3)*0.25*a3*(z[-1]**2)*(x[j]**3) + 0.5*const*(x[j]**2) + constx1*x[j]
        endx2 = (1/2)*a0*z[-1]*(x[j+1]**2) + (1/3)*a1*0.5*z[-1]*(x[j+1]**3) + (0.25*a2*z[-1]**2)*(x[j+1]**2) + (1/3)*0.25*a3*(z[-1]**2)*(x[j+1]**3) + 0.5*const*(x[j+1]**2) + constx1*x[j+1]   
        
        if j == 0:
            constx2 = 0 
        if j != 0:
            constx2 = xlist2[-1] - startx2
        startx2 += constx2
        endx2 += constx2
        
        xlist2.append(endx2)

        startx3 = (1/6)*a0*z[-1]*(x[j]**3) + (1/24)*a1*z[-1]*(x[j]**4) + ((1/12)*a2*z[-1]**2)*(x[j]**3) + (1/48)*a3*(z[-1]**2)*(x[j]**4) + (1/6)*const*(x[j]**3) + (1/2)*constx1*x[j]**2 + constx2*x[j]
        endx3 = (1/6)*a0*z[-1]*(x[j+1]**3) + a1*(1/24)*z[-1]*(x[j+1]**4) + ((1/6)*a2*z[-1]**2)*(x[j+1]**3) + (1/48)*a3*(z[-1]**2)*(x[j+1]**4) + (1/6)*const*(x[j+1]**3) + (1/2)*constx1*x[j+1]**2 + constx2*x[j+1]
        
        if j == 0:
            constx3 = 0
        if j != 0:
            constx3 = xlist3[-1] - startx3
        endx3 += constx3
        
        startx3 += constx3 
        xlist3.append(endx3) 
       
        sumert1 += endx3 - startx3
        sumert2 += endx1 - startx1
    print(sumert2)
    
            #Integrated function of x wrt x
    '''elif choice == 3:
            G = 0.5*part1*(x1**2) + (1/6)*part2*(x1**3) + 0.5*part3*(x1**2) + (1/6)*part4*(x1**3)
            integrsmall.append(G3)
            #Integrated function of x wrt x twice
        elif choice == 4:
            G = (1/6)*part1*(x1**3) + (1/24)*part2*(x1**4) + (1/6)*part3*(x1**3) + (1/24)*part4*(x1**4)
            integrsmall.append(G4)
            #Integrated function of x wrt x thrice
        elif choice == 5:
            G = (1/24)*part1*(x1**4) + (1/120)*part2*(x1**5) + (1/24)*part3*(x1**4) + (1/120)*part4*(x1**5)
            integrsmall.append(G5)
            #Integrated function of x wrt x four times
        else: 
            print("That is not a valid choice")
            print("Too bad, try again")
            break        '''
            
            
    return

integrFunc(1.611, coefarray, x, z, 1)

                

            
    
    
def cg(x2):
    intcg = []
    
    for j in range(len(x)-1):
        dothing = True
        integrcg = []
        for i in range(len(z)-1):
            a0 = float(coefarray[i,j,0])
            a1 = float(coefarray[i,j,1])
            a2 = float(0.5*coefarray[i,j,2])
            a3 = float(0.5*coefarray[i,j,3])
            
            part1 = a0*(z[i+1]-z[i])
            part2 = a1*(z[i+1]-z[i])
            part3 = a2*(z[i+1]**2-z[i]**2)
            part4 = a3*(z[i+1]**2-z[i]**2)
            
            CG = ((0.5*a0*(z[i+1]**2 - z[i]**2) + 0.5*a1*(z[i+1]**2 - z[i]**2)*x2 + 
                 (1/3)*a2*(z[i+1]**3-z[i]**3) + (1/3)*a3*(z[i+1]**3-z[i]**3)*x2)/(part1
                 + (x2)*part2 + part3 + (x2)*part4))
#            print(CG)
            
            if dothing == True:
                dothing = False
            
            integrcg.append(CG)
            
        integrcg = np.array(integrcg)
        integrcg = np.sum(integrcg)
            
        intcg.append(integrcg)
            
    intcg = np.array(intcg, dtype = 'f8')
         
    return intcg


f = interpolate.interp2d(x, z, q, kind = 'linear')
xnew = np.linspace(x[0], x[-1], 1000)
znew = np.linspace(z[0], z[-1], 1000)

qnew = f(xnew, znew)

print(trapz(trapz(qnew,znew), xnew))
#value = input('value of x: ')

value = x[40]  #          "DONT FORGET TO FIX"
print()

for i in range(len(x)-1):
    if float(value) > float(x[i]):       #determining in which part of the aileron we are and which g(x) to use
        pass
    else:
        t = g(value)
        print(t)
        finalg = t[i]
        print(finalg)                   #magintude of aerodynamic load on chord length for any x
#        break
    
        u = cg(value)
        finalcg = u[i]
        print('cg location resulatantn force is: ', finalcg)
        break


""" Creating test for the intergration function, by using array of ones"""


"""NUMERICAL VERIFICATION OF INTEGRATION:"""
        


