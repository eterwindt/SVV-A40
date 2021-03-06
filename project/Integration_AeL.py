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

#----------------------- functions definitions --------------------------------
Ca = 0.505 #[m], chord of aileron
la = 1.611 #[m], span of aileron

#-------------------------- definitions- --------------------------------------


#OBTAINING DATA
def readfile(filename):
    file = open(filename, "r")
    lines = file.readlines()                    # read all lines, returns list
    for n in range(0,len(lines)):
        lines[n] = lines[n].split(",")          # Split lines at "," to get individual entries
        lines[n][-1] = lines[n][-1][:-1]        # Remove \n at last entries
    file.close()
    
    return lines

#OBTAINING COORDINATES
def load_data(filename, Ca, la):
    q = readfile(filename)             #Importing all aerodynamic loading data into lists in lists
    Nx = len(q[0])
    Nz = len(q)
    
    zmin       = []
    x       = []
    for i in range(1,Nz+1):
        theta_z_i = (i-1)*pi/Nz
        theta_z_i1 = (i)*pi/Nz
    
        z_i = - 0.5*((Ca/2)*(1-cos(theta_z_i)) + (Ca/2)*(1-cos(theta_z_i1)))
        zmin.append(z_i)
    
    for j in range(1,Nx+1):
        theta_x_j = (j-1)*pi/Nx
        theta_x_j1 = (j)*pi/Nx
    
        x_j = 0.5*((la/2)*(1-cos(theta_x_j)) + (la/2)*(1-cos(theta_x_j1)))
        x.append(x_j)
    
    data = np.array(q, dtype='f8', copy=True, order='K')     #Coverting q into array
    
    test_data = np.ones_like(data)      #test case for the num model
    
    z = (-1)*np.array(zmin) #making data conform with macauley and section analysis
    data = np.fliplr(data)
#    data = 1000*data
    
    return data, test_data, x, z

#NUMERICAL INTERPOLATION
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

            intarray = [[1, float(x[i]), float(z[j]), float(x[i])*float(z[j])],
                        [1, float(x[i]), float(z[j+1]), float(x[i])*float(z[j+1])],
                        [1, float(x[i+1]), float(z[j]), float(x[i+1])*float(z[j])],
                        [1, float(x[i+1]), float(z[j+1]), float(x[i+1])*float(z[j+1])]]
                        # bivariate linear spline interpolation matrix
            dataarray = [[float(data[j,i])],  
                         [float(data[j+1,i])],
                         [float(data[j,i+1])],
                         [float(data[j+1,i+1])]]
                        # solution matrix from aerodynnamic loading
            coef = np.dot(np.linalg.inv(intarray), dataarray)  
            # obtaining coefficient matrix for solving linear splines
            
            coefsmall.append(coef)
            # creating a matrix with coefficients for all the data points on the aileron in x-dir
            if dothing == True:
                # print(coef)
                dothing = False

        coefsmall = np.array(coefsmall)
        coefsmall = np.reshape(coefsmall, (40, 1, 4, 1))

        coefarray.append(coefsmall)
    coefarray = np.array(coefarray, dtype='f8')  
    # creating the matrix in z-dir by appending the x-dir matrices
    
    coefarray = np.reshape(coefarray, (80, 40, 4))
    return coefarray

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

#NUMERICAL INTEGRATION
def print_menu():
    print( "Choose an option")
    print() 
    print( "1) Find Rforce at an x-coordinate due to aerodynamic loading")
    print( "2) Force f(x) integrated once, up to x")
    print( "3) Force f(x) integrated twice, up to x")
    print( "4) Force f(x) integrated thrice, up to x")
    print( "5) Force f(x) integrated four times, up to x")


def g(x1, coefarray, choice, x, z):
    '''
    :param x1: location at aileron at which or to where to integrate
    :coefarray: The coefficients obtained for interpolation
    :choice: How many times to integrate (see menu for options)
    :param x: list of x coords
    :param z: list of z coords
    :return: value integration, resultant force per square
    '''
    """Integration"""
    integr = []
    for j in range(len(x)-1):
        dothing = True
        integrsmall = []
        integrcgz = []
        for i in range(len(z)-1):
            a0 = float(coefarray[i,j,0])        
            a1 = float(coefarray[i,j,1])
            a2 = float(0.5*coefarray[i,j,2])
            a3 = float(0.5*coefarray[i,j,3])
            
            part1 = a0*(z[i+1]-z[i])
            part2 = a1*(z[i+1]-z[i])
            part3 = a2*(z[i+1]**2-z[i]**2)
            part4 = a3*(z[i+1]**2-z[i]**2)
           
            
            cgz = part1 + float(x1)*part2 + part3 + float(x1)*part4
            integrcgz.append(cgz)
            
            if choice == 1:
                G = part1 + float(x1)*part2 + part3 + float(x1)*part4
                integrsmall.append(G)
                #Integrated function of aerodynamic loading, wrt z to get f(x)               
            elif choice == 2:
                if float(x1)>=x[j+1]:
                    G = part1*(x[j+1]-x[j]) + 0.5*part2*(x[j+1]**2-x[j]**2) +\
                    part3*(x[j+1]-x[j]) + 0.5*part4*(x[j+1]**2-x[j]**2)
                    integrsmall.append(G)
                else:
                    G = part1*(float(x1)-x[j]) + 0.5*part2*(float(x1)**2-x[j]**2) +\
                    part3*(float(x1)-x[j]) + 0.5*part4*(float(x1)**2-x[j]**2)
                    integrsmall.append(G)
                    print("test works")
                    break
                #Integrated function of x wrt x
            elif choice == 3:
                if float(x1)>=x[j+1]:
                    G = 0.5*part1*(x[j+1]**2-x[j]**2) + (1/6)*part2*(x[j+1]**3-x[j]**3) +\
                    0.5*part3*(x[j+1]**2-x[j]**2) + (1/6)*part4*(x[j+1]**3-x[j]**3)
                    integrsmall.append(G)
                else:
                    G = 0.5*part1*(float(x1)**2-x[j]**2) + (1/6)*part2*(float(x1)**3-x[j]**3) +\
                    0.5*part3*(float(x1)**2-x[j]**2) + (1/6)*part4*(float(x1)**3-x[j]**3)
                    integrsmall.append(G)
                    break
                #Integrated function of x wrt x twice
            elif choice == 4:
                if float(x1)>=x[j+1]:
                    G = (1/6)*part1*(x[j+1]**3-x[j]**3) + (1/24)*part2*(x[j+1]**4-x[j]**4) +\
                    (1/6)*part3*(x[j+1]**3-x[j]**3) + (1/24)*part4*(x[j+1]**4-x[j]**4)
                    integrsmall.append(G)
                else:
                    G = (1/6)*part1*(float(x1)**3-x[j]**3) + (1/24)*part2*(float(x1)**4-x[j]**4) +\
                    (1/6)*part3*(float(x1)**3-x[j]**3) + (1/24)*part4*(float(x1)**4-x[j]**4)
                    integrsmall.append(G)
                    break
                #Integrated function of x wrt x thrice
            elif choice == 5:
                if float(x1)>=x[j+1]:
                    G = (1/24)*part1*(x[j+1]**4-x[j]**4) + (1/120)*part2*(x[j+1]**5-x[j]**5) +\
                    (1/24)*part3*(x[j+1]**4-x[j]**4) + (1/120)*part4*(x[j+1]**5-x[j]**5)
                    integrsmall.append(G)
                else:
                    G = (1/24)*part1*(float(x1)**4-x[j]**4) + (1/120)*part2*(float(x1)**5-x[j]**5)+\
                    (1/24)*part3*(float(x1)**4-x[j]**4) + (1/120)*part4*(float(x1)**5-x[j]**5)
                    integrsmall.append(G)
                    break
                #Integrated function of x wrt x four times
            else: 
                print("That is not a valid choice")
                print("Too bad, try again")
                break
            if dothing == True:
                dothing = False
        integrsmall = np.array(integrsmall)
        integrcgz = np.array(integrcgz)
        
        resultants_z_per_square = integrcgz
        
        integrsmall = np.sum(integrsmall) 
        #summing the parts in the z direction to gain the force in x
        
        integr.append(integrsmall)
    integr = np.array(integr, dtype = 'f8') #array containing integrated functions to get g(x)
    return integr, resultants_z_per_square


#%% Test cases: 
    

filename = 'aerodynamicloadf100.dat'    

data, test_data, x, z = load_data(filename, Ca, la)



""" verifying the integration
veritestx = [0,1,2,3,4,5]
veritestz = [0,0.5,1.0,1.5,2.0]
vericoef = np.ones((5,6,4))
testchoice = 2

verification = g(5, vericoef, testchoice, veritestx, veritestz)
answer = verification[0]
if testchoice==1:
    print("expected answers:")
    print("x = 1 -> 8.0")
    print("x = 3 -> 16")
    print("x = 5 -> 24")
    print("test = ", answer[0])
if testchoice>1:
    if testchoice==2:
        print("expected answer for full range integral: 70")
    if testchoice==3:
        print("expected answer for full range integral: 133.33")
    if testchoice==4:
        print("expected answer for full range integral: 187.5")
    if testchoice==5:
        print("expected answer for full range integral: 208.33")
    print("test = ", sum(answer))
    

The test integration the same as the handcalculated values so the integration model is correct.
"""


"""
Verifying the z-loaction (choice==1) and x-location (choice>1) of the resultant force

testmode = input("Test mode? (y/n): ")      #entering test case to verify num model

if testmode.lower() == 'y':
    data = test_data

From this we see that our Fres location calculation is correct, 
since for a linear distributed load,the Fres acts in the center of the aileron. 
This is also the outcome of our model with the test case"""

#%% Main Program:

#print_menu()
#choice = int(input("Your choice: "))
#print()

#value = input('value of x: ')
value = 1.608
choice = 2      #set to one if you want the location of Fres
coeffs = coefficients(data, x, z)

"""for the Fres per slice and location of Fres"""
#if choice ==1:    
finallist = []
for i in range(len(x)-1):
    if float(value) > float(x[i]):       
        pass
        #determining in which part of the aileron we are and which g(x) to use
        
    else:
        t, resultants_z_per_square = g(value, coeffs, choice, x, z)
        finalg = t[i]
        finallist.append(finalg)
        if choice ==1:
            print("Final integrated value = ", finalg)  
        save = i
        break
        #magintude of aerodynamic load on chord length for any x integrated per your choice

cglistz = []                       #getting the resultant force arm by a = F/M
cglistx = []
for i in range(len(z)-1):
    centerz = (z[i+1] + z[i])/2
    cglistz.append(centerz)
    np.array(cglistz)
for j in range(len(x)-1):
    centerx = (x[j+1] + x[j])/2
    cglistx.append(centerx)
    np.array(cglistx)
    
"""For the interation of forces """ 
#if choice >1:
final = g(value, coeffs, choice, x, z)   #edit 3rd input to change number of integrations
print("Integrated value: ", np.sum(final[0]))

"""Locations of Fres"""
   
if choice ==1:
    zloc_resultantforce_perx = sum(cglistz*resultants_z_per_square)/sum(resultants_z_per_square)
    print("z location resultant force is: ", zloc_resultantforce_perx)
if choice >1:
    xloc_resultantforce = sum(cglistx*final[0])/sum(final[0])
    print("x location Fres full aileron is:", xloc_resultantforce) #for full aileron set x>=1.611




