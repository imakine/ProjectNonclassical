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

    # Initial parameters

    initial = time.time()	 # Initial time
    sigt1=1.0			 # Cross section of the solid
    l1 = 1.0;			 # Length of the solid media
    l2 = 0;			 # Length of the void media
    l = l1+l2;			 # Total period of the media
    n=10000;  			 # Number of particles
    N= 16 			 # Number of direction



    # Computation of the different value of mu (cosine of theta) for each direction and its weight in the Gauss-Legendre approximation

    mu_list = np.polynomial.legendre.leggauss(N)[0][int(N/2):N] #use only positive directions
    weight_list = np.polynomial.legendre.leggauss(N)[1][int(N/2):N] #use only positive directions

    # Initialisation of the final result for the different moment of s (until 6th)

    moys1=moys2=moys3=moys4=moys5=moys6 = 0;

    # Loop on each direction of the function giving the distance travelled by the particle
    for i in range(int(N/2)):

	# Print some intermediary results, number of the iteration and the value of mu

        print("Direction nr : ",i+1)
        print("Direction mu = ",mu_list[i])
        dist1,dist2,dist3,dist4,dist5,dist6 = num_mc(n,mu_list[i],l1,l2,sigt1);

	# Computation of the first moment, print of its value and add with the weight to the final result. Computation of the standart deviation

        mm1=np.mean(dist1)
        print("moy s1 for this direction = ", mm1)
        st1=np.std(dist1) #This is the standard deviation
#        st1=(np.std(dist1))**2+(np.mean(dist1))**2;  #I don't know what you are doing in this line
        moys1 = moys1 + weight_list[i]*mm1

	# Computation of the second moment, print of its value and add with the weight to the final result. Computation of the standart deviation

        mm2=np.mean(dist2)
        print("moy s2 for this direction = ", mm2)
        st2=np.std(dist2)
        #st2=(np.std(dist2))**2+(np.mean(dist2))**2;
        moys2 = moys2 + weight_list[i]*mm2

	# Computation of the third moment, print of its value and add with the weight to the final result. Computation of the standart deviation

        mm3=np.mean(dist3)
        print("moy s3 for this direction = ", mm3)
        st3=np.std(dist3)
        #st3=(np.std(dist3))**2+(np.mean(dist3))**2;
        moys3 = moys3 + weight_list[i]*mm3

	# Computation of the fourth moment, print of its value and add with the weight to the final result. Computation of the standart deviation

        mm4=np.mean(dist4)
        print("moy s4 for this direction = ", mm4)
        st4=np.std(dist4)
        #st4=(np.std(dist4))**2+(np.mean(dist4))**2;
        moys4 = moys4 + weight_list[i]*mm4

	# Computation of the fifth moment, print of its value and add with the weight to the final result. Computation of the standart deviation

        mm5=np.mean(dist5)
        print("moy s5 for this direction = ", mm5)
        st5=np.std(dist5)
        #st5=(np.std(dist5))**2+(np.mean(dist5))**2;
        moys5 = moys5 + weight_list[i]*mm5

	# Computation of the sixth moment, print of its value and add with the weight to the final result. Computation of the standart deviation

        mm6=np.mean(dist6)
        print("moy s6 for this direction = ", mm6)
        st6=np.std(dist6)
#        st6=(np.std(dist6))**2+(np.mean(dist6))**2;
        moys6 = moys6 + weight_list[i]*mm6

    final = time.time()

    # Print of the final results for the 6 first moments

    print("Computation of the moments of 's' for " + str(N) + " directions")
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

# This function gives the distance travelled by a particle for a given direction
# n is the number of particles
# mu is the considered direction
# l1 is the length of the first solid media
# l2 is the length of the second void media
# sigmq1 is the cross section in the first media


def num_mc(n,mu,l1,l2,sigma1):

    # Initialisation

    m=0;
    dist=np.zeros(n);
    dist2=np.zeros(n);
    dist3=np.zeros(n);
    dist4=np.zeros(n);
    dist5=np.zeros(n);
    dist6=np.zeros(n);

    while m<=n-1:

	# First random generation of the initial position in the solid media

        x0=l1*random.random();

	# Sampling of the exponential

        s1=-math.log(1-random.random())/sigma1;

        if s1<= x0/mu:	# This case represents if the particle stays in the first solid media or go to the next
            s=s1;
        else:		# The particle goes to the next solid media (or more far than that)
            s=s1+(l2/mu)*math.ceil(mu*(s1-x0/mu)/l1);
        dist[m]=s;
        dist2[m]=s**2;
        dist3[m]=s**3;
        dist4[m]=s**4;
        dist5[m]=s**5;
        dist6[m]=s**6;
        m=m+1;
    return list(dist),list(dist2),list(dist3),list(dist4),list(dist5),list(dist6)

main()
