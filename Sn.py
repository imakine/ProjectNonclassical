# -*- coding: utf-8 -*-
"""
Created on Thu May 18 01:28:26 2017

@author: Ilker
"""
import numpy as np
import scipy.optimize
import random
import math
import time

def main():
    initial = time.time()
    sigt1=1.0 # cross section
    sigt2 = 0.0
    l1 = 1.0;
    l2 = 1.0;
    l = l1+l2;
    n=100000;  # number of particles
    N= 8   
    
    mu_list = np.polynomial.legendre.leggauss(N)[0]
    weight_list = np.polynomial.legendre.leggauss(N)[1]
    
    moys1=moys2=moys3=moys4=moys5=moys6 = 0;
    for i in range(N):
        print("Direction nr : ",i+1)
        print("Direction mu = ",mu_list[i])
        dist1,dist2,dist3,dist4,dist5,dist6 = num_mc(n,mu_list[i],l1,l2,sigt1); 
                      
        mm1=np.mean(dist1)
        print("moy s1 for this direction = ", mm1)
        st1=(np.std(dist1))**2+(np.mean(dist1))**2;
        moys1 = moys1 + weight_list[i]*mm1  
                                   
        mm2=np.mean(dist2)
        print("moy s2 for this direction = ", mm2)
        st2=(np.std(dist2))**2+(np.mean(dist2))**2;
        moys2 = moys2 + weight_list[i]*mm2
                                   
        mm3=np.mean(dist3)
        print("moy s3 for this direction = ", mm3)
        st3=(np.std(dist3))**2+(np.mean(dist3))**2;
        moys3 = moys3 + weight_list[i]*mm3
                                   
        mm4=np.mean(dist4)
        print("moy s4 for this direction = ", mm4)
        st4=(np.std(dist4))**2+(np.mean(dist4))**2;
        moys4 = moys4 + weight_list[i]*mm4
                                   
        mm5=np.mean(dist5)
        print("moy s5 for this direction = ", mm5)
        st5=(np.std(dist5))**2+(np.mean(dist5))**2;
        moys5 = moys5 + weight_list[i]*mm5
        
        mm6=np.mean(dist6)
        print("moy s6 for this direction = ", mm6)
        st6=(np.std(dist6))**2+(np.mean(dist6))**2;
        moys6 = moys6 + weight_list[i]*mm6
                                   
    final = time.time() 

    print("Computation of the moments of 's' for " + str(2*N) + " directions")
    print(' ')                              
    print("First moment = " + str (moys1) + " with a st = " + str(st1));
    print(" "); 
    print("2nd moment = " + str (moys2) + " with a st = " + str(st2));
    print(" ");
    print("3nd moment = " + str (moys3) + " with a st = " + str(st3));
    print(" ");     
    print("4th moment = " + str (moys4) + " with a st = " + str(st4));
    print(" ");
    print("5th moment = " + str (moys5) + " with a st = " + str(st5));
    print(" ");
    print("6th moment = " + str (moys6) + " with a st = " + str(st6));
    print(" ");   
    print(" Time elapsed during the running of the code : ", final - initial)     
def num_mc(n,mu,l1,l2,sigma1):
    m=0; 
    dist=np.zeros(n);
    dist2=np.zeros(n);
    dist3=np.zeros(n);
    dist4=np.zeros(n);
    dist5=np.zeros(n);
    dist6=np.zeros(n);
    while m<=n-1:
        x0=l1*random.random();
        s1=-math.log(1-random.random())/sigma1; 
        if s1<= x0/mu:
            s=s1;
        else:
            s=s1+(l2/mu)*math.ceil(mu*(s1-x0/mu)/l1);
        dist[m]=s;
        dist2[m]=s**2;
        dist3[m]=s**3;
        dist4[m]=s**4;    
        dist5[m]=s**5;
        dist6[m]=s**6;    
        m=m+1;
    return list(dist),list(dist2),list(dist3),list(dist4),list(dist5),list(dist6)    
