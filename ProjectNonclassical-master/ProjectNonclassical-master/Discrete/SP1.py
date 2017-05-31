# -*- coding: utf-8 -*-
"""
Created on Thu May  4 01:04:26 2017

@author: Ilker
"""
from numpy import *
import matplotlib.pyplot as plt
from numpy import linalg as LA
import numpy as np
import math

def main():
    
    # Input parameters
    
    T = int(input("Total length of the system : "));
    n = int(input("Number of points : "));
    m1 = float(input("First moment of s : "));
    m2 = float(input("Second moment of s : "));
    c = float(input("Scattering ratio c : "))        
    q = float(input("Source : "));                
    
    #Discretization & calculation

    h = T/n; #step
    sigmatot = 2*m1/m2;
    sigmaabs = (1-c)/m1;
    D = 1/(3*sigmatot);

    #initialisation of matrices
    
    A = np.zeros((n+1,n+1));
    b = np.zeros((1,n+1));
    x = np.zeros((1,n+1));         
        
    # Coefficients and constants

    A[0,0] = 1/2 + 3*D/(2*h);
    A[0,1] = -4*D/(2*h);
    A[0,2]= D/(2*h);
    print(A)
    for i in range(1,n):
        A[i,i-1] = -D/h**2;
        A[i,i] = sigmaabs + 2*D/h**2;
        A[i,i+1] = -D/h**2;
        b[0,i] = q;
    print(A)  
    print("bb = ", b)
    A[n,n-2]=D/(2*h);
    A[n,n-1]=-4*D/(2*h);
    A[n,n]=1/2 + 3*D/(2*h);
    print(A)
    solution = gauss(A,b,x,2000);
    print(solution.tolist())
    """
    h1 = int(1/h);
    print(h1)
    xaxis = [0]*h1;
    print(xaxis)
    sumstep = 0;        
    while sumstep <len(xaxis):
        xaxis[i]=h; 
        sumstep = sumstep + h; 
    """             
    i=-T/2;
    j = 0;
    xaxis = [];
    xref = []        
    while i < T/2+h: 
        xref.append(i);
        xaxis.append(j)           
        i = i + h;
        j = j + h;
        
    reference = [0]*len(xaxis);    
    L = math.sqrt(D/sigmaabs)  
    for i in range(len(reference)):
        reference[i] = (q/sigmaabs)*(1- math.cosh(xref[i]/L)/(math.cosh(T/L)+(D/L)*math.sinh(T/L)))  
    print(reference)    
    plt.plot(xaxis,list(solution[0]));
    plt.plot(xaxis,list(reference));
    plt.ylabel('Flux(z)');
    plt.show(); 

def gauss(A, b, x0, maxiter):
    tol = 1e-6
    k = 1
    n = len(A)
    x = x0;
    if n != len(A[0]):
        print("Matrix must be n x n")
    while k < maxiter:
        for i in range(n):
            sums = 0;
            for j in range(n):
                if i != j:
                    sums = sums + A[i,j]*x[0,j]
            x[0,i]= (b[0,i]-sums)/A[i,i]
        if ( LA.norm(x0-x) < tol):
                return x
        k = k +1;
        for i in range(n):
            x0[0,i] = x[0,i]
    return x        
            
            
    