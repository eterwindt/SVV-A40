# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 13:26:35 2020

@author: esmee
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:56:53 2020

@author: esmee
"""
#----------------------------- imports ----------------------------------------
from numpy import pi, cos
import numpy as np
from matplotlib import pyplot as plt

#----------------------- functions definitions --------------------------------
Ca = 0.505 #[m], chord of aileron
la = 1.611 #[m], span of aileron

#-------------------------- main program --------------------------------------



def readfile(filename):
    file = open(filename, "r")
    lines = file.readlines()                    # read all lines, returns list
    for n in range(0,len(lines)):
        lines[n] = lines[n].split(",")          # Split lines at "," to get individual entries
        lines[n][-1] = lines[n][-1][:-1]        # Remove \n at last entries
    file.close()
    
    return lines

def load_data(filename):
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
    
    data = np.array(q, dtype='f8', copy=True, order='K')     #Coverting q into array
    
    test_data = np.ones_like(data)      #test case for the num model
    
    return data, test_data, x, z


def coefficients(data, x, z):
    '''
    :param data: loading data
    :param x: list of x coords
    :param z: list of z coords
    :return: coefficients
    '''
    """ Now interpolate wrt to x and y"""
    coefarray = []
    for j in range(len(z) - 1):
        dothing = True
        coefsmall = []
        for i in range(len(x) - 1):

            intarray = [[1, float(x[i]), float(z[j]), float(x[i]) * float(z[j])],
                        # bivariate linear spline interpolation matrix
                        [1, float(x[i]), float(z[j + 1]), float(x[i]) * float(z[j + 1])],
                        [1, float(x[i + 1]), float(z[j]), float(x[i + 1]) * float(z[j])],
                        [1, float(x[i + 1]), float(z[j + 1]), float(x[i + 1]) * float(z[j + 1])]]

            dataarray = [[float(data[j, i])],  # solution matrix from aerodynnamic loading
                         [float(data[j + 1, i])],
                         [float(data[j, i + 1])],
                         [float(data[j + 1, i + 1])]]
            coef = np.dot(np.linalg.inv(intarray), dataarray)  # obtaining coefficient matrix for solving linear splines

            coefsmall.append(
                coef)  # creating a matrix with coefficients for all the data points on the aileron in x-dir
            if dothing == True:
                # print(coef)
                dothing = False

        coefsmall = np.array(coefsmall)
        coefsmall = np.reshape(coefsmall, (40, 1, 4, 1))

        coefarray.append(coefsmall)
    coefarray = np.array(coefarray, dtype='f8')  # creating the matrix in z-dir by appending the x-dir matrices
    coefarray = np.reshape(coefarray, (80, 40, 4))
    return coefarray
    # print(coefarray[0, 0, 0:4])


#coefarray = coefficients(test_data, x, z)

"""
Numerical Verification of the linear bivariate spline interpolation by  
comparing the outcome of the function with the original data           
for j in range(len(z)-1):
        for i in range(len(x)-1):
            coefs = coefarray[j,i,:]           
            function = coefs[0] + coefs[1]*x[i] + coefs[2]*z[j] + coefs[3]*x[i]*z[j]
#            print (function)
This holds true, since the output of the test is conform with the input
"""     
        
"""Now to integrate the function wrt z to obtain a 1D f(x),
and integrate wrt x"""

def print_menu():
    print( "Choose an option")
    print() 
    print( "1) Find Rforce at an x-coordinate due to aerodynamic loading")
    print( "2) Force f(x) integrated once")
    print( "3) Force f(x) integrated twice")
    print( "4) Force f(x) integrated thrice")
    print( "5) Force f(x) integrated four times")


def g(x1, coefarray, choice):
    integr = []
    
    for j in range(len(x)-1):
        dothing = True
        integrsmall = []
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
           
            
            cg = part1 + (x1)*part2 + part3 + (x1)*part4
            integrcg.append(cg)
            
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
        integrcg = np.array(integrcg)
        resultants_z_per_square = integrcg
        integrsmall = np.sum(integrsmall) #summing the parts in the z direction to gain the force for a coordinate x
        
        integr.append(integrsmall)
       
    integr = np.array(integr, dtype = 'f8') #array containing integrated functions to get g(x)
    
    return integr, resultants_z_per_square



#%% Main Code:    

filename = 'aerodynamicloadf100.dat'    

data, test_data, x, z = load_data(filename)

testmode = input("Test mode? (y/n): ")      #entering test case to verify num model

if testmode.lower() == 'y':
    data = test_data

print_menu()
choice = int(input("Your choice: "))
print()
"""From this we see that our Fres location calculation is correct, since for a linear distributed load,
the Fres acts in the center of the aileron. This is also the outcome of our model with the test case"""

#value = input('value of x: ')
value = x[0]  #          "DONT FORGET TO FIX"
print()

coeffs = coefficients(data, x, z)

for i in range(len(x)-1):
    if float(value) > float(x[i]):       #determining in which part of the aileron we are and which g(x) to use
        pass
    else:
        t, resultants_z_per_square = g(value, coeffs, choice)
        finalg = t[i]
        print("Final integrated value = ", finalg)  #magintude of aerodynamic load on chord length for any x
        save = i                                    #integrated per your choice
        break
    
cglist = []                       #getting the resultant force arm by a = F/M
for i in range(len(z)-1):
    center = (z[i+1] + z[i])/2
    cglist.append(center)
    np.array(cglist)
   
zloc_resultantforce_perx = sum(cglist*resultants_z_per_square)/sum(resultants_z_per_square)
print("cg location resultant force is: ", zloc_resultantforce_perx)

#%%
"""Testing if integration choice==1 is correct """
from scipy.integrate import trapz

def calc_error(real, approx):           #relative error
    error = real - abs(approx)
    return np.abs(error) / real

datatest = data[:,save]

intgration_ref = trapz(z, datatest)
print("The integration test = ",intgration_ref)
error = calc_error(intgration_ref, finalg)
print('The error = ', error)








