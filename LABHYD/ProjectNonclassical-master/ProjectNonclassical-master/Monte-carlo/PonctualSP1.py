# -*- coding: utf-8 -*-
"""
Created on Mon May 15 21:19:50 2017

@author: Ilker
"""

import random
import math
import time
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize 

def main(c):
    ############################################################
    # INPUT
    #
    # c : the scattering ratio of the media
    #
    #
    # OUTPUT
    # Flux SP123
    #
    #############################################################

    #Initialisation of the parameters of the simulation
    thickness = 100;                    #Thickness of the mediq
    sigmaabs = 1-c;                     #Absorption cross-section of the media
    sigmascat = c;                      #Scattering cross-section of the media
    sigmatot=sigmaabs+sigmascat;        #Total cross-section of the media
    global lamb;                        #Lamb is the mean free path
    if (sigmatot !=0):
        lamb = 2/((3)**(0.5)*sigmatot); #mean-free path of SP1
    Q = 1                               #Source of neutrons
    n = 100000;                          #Number of path
    step = 0.00275;                       #Step for bins Number of bins will be Length/step

    #Output file
    filename = 'realend';
    file = open(filename + ".txt", "w") 
    fluxPSP1 = open("fluxPSP1.txt","w")
    coordPSP1 = open("coordPSP1.txt","w")
    variancePSP1 = open("variancePSP1.txt","w");
    deviationPSP1 = open("deviationPSP1.txt","w");
    
    # Time starts

    initialtime = time.time();
   

    # All parameters to count
    
    global numlost;                                                     #Number of particles lost by the left side (z < 0)
    global numesc;                                                      #Number of particle lost by the rifht side (z > Thickness)
    global numdeath;                                                    #Number of particles absorbed in the media
    global numscattot;                                                  #Number of collisions for ALL particles
    global freepath;                                                    #The distance giving by the sampling
    global flux_local;                                                  #List of an estimator of the flux in the media ; collision between 0 and 1 is in first place, between 1 and 2 is second ...
    global s;                                                           #List of the 6 first moments of s                                                
    global variancepic;                                                 #List of the sum of the square of each scattering                              
    global deviationpic;                                                #List of the sum of each scattering
    test = [0]*(int(thickness/step));                                   #List of the variance of each element in flux_local 
    test2 = [0]*(int(thickness/step));                                  #List of the standart deviation of each element in flux_local
    fluxpic = [0]*(int(thickness/step));
    variancepic = [0]*(int(thickness/step));
    s =[0]*6;
    flux_local = [0]*(int(thickness/step));  
    numlost=0;
    numesc=0;
    numdeath =0;
    numscattot = 0;
    freepath =0;
    
    # Loop on each history
    
    for i in range(1,n+1):
        # Interface with user
        if (i%1000 == 0):
            print(str(i/n*100) + "% of the running code done")
        randnum = random.random();
        z0=thickness/2+(2*randnum-1)/2;                                 #Initial rqndom position of particles according to the regional source : Source = Q between -0.5 + thickness/2 and +0.5 + thickness/2
        randnum = random.random()
        theta = 180*randnum;                                            #Initial random angle for the particle
        w0=math.cos(math.pi/180*theta);                                 #Initial direction of the particle
        z=z0;                                                           # position z follows the particle
        w=w0;                                                           # direction w follows the particle 
        ilost=0;                                                        #Event escape by left side value 0 --> particle does not escape by the left side, value 1 --> particle has escaped by the left side
        iesc=0;                                                         #Event escape by right side value 0 --> particle does not escape by the right side, value 1 --> particle has escaped by the right side                                    
        ideath=0;                                                       #Event absorption by the media value 0 --> particle is not absorpbed, value 1 --> particle has been absorbed 
        iscat=0;                                                        #Number of collision by particle                
        ivariancepic = [0]*(int(thickness/step));                       #Number of collision squared by particle
        iflux = [0]*(int(thickness/step));                              #Number of collision by particle 
        
        #loop on the life of the particle 
        
        while((ilost+iesc+ideath)==0):
            randnum=random.random();               
            if (sigmatot !=0):
                freepath = scipy.optimize.brentq(fsp,0.0,100.0,args=(randnum))/((3**0.5)*sigmatot) #Sampling of the distance travelled
                s[0] = s[0] + freepath;                                 #First moment estimator
                s[1] = s[1] + freepath**2;                              #Second moment estimator
                s[2] = s[2] + freepath**3;                              #Third moment estimator    
                s[3] = s[3] + freepath**4;                              #..
                s[4] = s[4] + freepath**5;                              #..
                s[5] = s[5] + freepath**6;                              #.
            z1=z+w*freepath;                                            #Update of the position
            if(z1 >= thickness or sigmatot==0):                         #Check if the particle is lost in left side
                iesc=1;
            elif(z1 < 0.):                                              ##Check if the particle is lost in right side                                              
                ilost=1; 
            else:
                xi = random.random();
                if (xi < sigmascat/sigmatot):                           #Check if it's a scattering with a probability of Sigmascat/sigmatot
                    index = int(z1/step)+1;                             #Count in the bins
                    flux_local[index-1] = flux_local[index-1] + 1;      #Increment event for flux
                    iflux[index-1] = iflux[index-1] +1;                 #Increment event for flux
                    ivariancepic[index-1] = ivariancepic[index-1] + 1;  #Increment event for flux        
                    iscat = iscat + 1;                                  #Increment of the number of scattering for the particle 
                    z=z1;                                               #Save the new position
                    randnum=random.random()                             #Isotropic collision, direction random
                    w=2*randnum-1;
                if (sigmaabs !=0 and xi > sigmascat/sigmatot):          #Check if it's an absorption with a probability Sigmaabs/sigmatot
                    index = int(z1/step) + 1;                           #Count in the bins
                    flux_local[index-1] = flux_local[index-1] + 1;      #Increment event for flux
                    iflux[index-1] = iflux[index-1] +1;                 #Increment event for flux
                    ivariancepic[index-1] = ivariancepic[index-1] + 1;  #Increment event for flux          
                    ideath = 1;                                         #Kill the particle
                else:
                    ideath = 0;                                         #No absorption 
        #Computation of the flux and the variance            
        for i in range(len(fluxpic)):            
            fluxpic[i] = fluxpic[i] + iflux[i];                         #We add the collision of each particle for flux          
            variancepic[i] = variancepic[i] + ivariancepic[i]**2;       #We add the square power of each particle for variance                                          
        numlost=numlost+ilost;                                          #Update of the number of lost in right side
        numesc=numesc+iesc;                                             #Update of the number of lost in left side
        numdeath=numdeath+ideath;                                       #Update of the number of qbsorption
        numscattot = numscattot + iscat;                                #Update of the number of scattering            
    for i in range(len(fluxpic)):
        variancepic[i] = variancepic[i]/(n-1);                          #Average number of collision on the sample
        fluxpic[i] = fluxpic[i]/n;                                      #Average number of collision square on the sample            
        test[i] = variancepic[i]-(fluxpic[i])**2;                       #Variance
        test2[i] = test[i]**(0.5);                                      #Standart deviation
    flux_local[:] = [(x/n)*(Q)/step/sigmatot for x in flux_local];      #Average flux with good normalisation The integral of the source appears Q*1 (because Q is constant in a distance = 1)
    s[:] = [x / (numscattot+numdeath+numesc+numlost) for x in s];       #Average of the moments     
    finaltime = time.time();
    print("n =",n);
    print("numscattot = ", numscattot);
    print("numlost =",numlost);
    print("numabs =",numdeath);
    print("numesc=",numesc);
    print("flux(z) = ",flux_local);
    print("s = ",s);   
    print("Time elapsed during the running of the code : ",finaltime - initialtime, "seconds");
    print(" ");
    print("Variance of flux = ",test)  
    print(" ");
    print('max of the flux = ', max(flux_local))
    print("average on bins max = ", (flux_local[flux_local.index(max(flux_local))]+flux_local[flux_local.index(max(flux_local))+1])/2)
    #print("Standart deviation of flux = ", test2)
    print(" ");    
    #print("Maximum for the standart deviation = "+ str(max(test2)) + " at the position " + str(test2.index(max(test2))) )
    plt.plot(flux_local)
    plt.ylabel('Flux(z)')
    plt.show()
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   ------------------------------------------------------------\n")
    file.write("   |                                                          |\n")
    file.write("   |                 MONTE-CARLO CODE FOR                     |\n")
    file.write("   |           SP1 EQUATION (1D-Ponctual source)              |\n")
    file.write("   |             Done at Berkeley by Ilker Makine             |\n")
    file.write("   |                      OUTPUT FILE                         |\n")
    file.write("   |                                                          |\n")
    file.write("   |                                                          |\n")
    file.write("   ------------------------------------------------------------\n")
    file.write("   ---------------------------INPUT PARAMETERS---------------------\n")  
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   " + " Thickness of the media = " + str(thickness) + "\n")
    file.write("   " + " Absorption cross-section = " + str(sigmaabs)+ "\n")
    file.write("   " + " Scattering cross-section = " + str(sigmascat)+ "\n")
    file.write("   " + " Total cross-section = " + str(sigmatot)+ "\n")
    file.write("   " + " Source of neutrons = " + str(Q)+ "\n")
    file.write("   " + " Number of particles = " + str(n)+ "\n")
    file.write("   " + " Precision in the media (step) = " + str(step)+ "\n")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   ---------------------------TIME ELAPSED---------------------\n")  
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   " + str(finaltime - initialtime) + " seconds are elapsed during the running of the code")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   ---------------------------MOMENTS--------------------------\n")  
    file.write("\n")
    file.write("\n")
    file.write("\n")
    s1 = 1.1547/sigmatot;
    s2 = 2/sigmatot**2;
    s3 = 5.6188/sigmatot**3;
    s4 = 13.3333/sigmatot**4;
    s5 = 48.188/sigmatot**5;
    s6 = 186.66667/sigmatot**6;
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^m         | Theoretical values | Computational values |        Errors        |\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^1         | "+str(s1)+" "*(22-len(str(s1))) +         " |"      +    str(s[0])+" "*(22-len(str(s[0])))+"|"          +str((s[0]- s1)/s1) + " "*(22-len(str((s[0]- s1)/s1)))+     "|\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^2         | "+str(s2)+" "*(22-len(str(s2))) +         " |"      +    str(s[1])+" "*(22-len(str(s[1])))+"|"        +str((s[1]- s2)/s2) +" "*(22-len(str((s[1]- s2)/s2))) +    "|\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^3         | "+str(s3)+" "*(22-len(str(s3))) +         " |"      +    str(s[2])+" "*(22-len(str(s[2])))+"|"           +str((s[2]- s3)/s3) +  " "*(22-len(str((s[2]- s3)/s3))) +      "|\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^4         | "+str(s4)+" "*(22-len(str(s4))) +         " |"      +    str(s[3])+" "*(22-len(str(s[3])))+"|"         +str((s[3]- s4)/s4)+    " "*(22-len(str((s[3]- s4)/s4)))  +     "|\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^5         | "+str(s5)+" "*(22-len(str(s5))) +         " |"      +    str(s[4])+" "*(22-len(str(s[4])))+"|"         +str((s[4]- s5)/s5)+  " "*(22-len(str((s[4]- s5)/s5)))   +      "|\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   |      s^6         | "+str(s5)+" "*(22-len(str(s5))) +         " |"      +    str(s[5])+" "*(22-len(str(s[5])))+"|"         +str((s[5]- s6)/s6)+ " "*(22-len(str((s[5]- s6)/s6))) +        "|\n")
    file.write("   |                  |                    |                      |                      |\n")
    file.write("   ---------------------------------------------------------------------------------------\n")
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   ---------------------------EVENTS COUNTER---------------------\n")  
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   " + str(numscattot) + " scatterings occurs during the simulation\n")
    file.write("   " + str(numlost+numesc) + " escapes occurs during the simulation\n")
    file.write("   " + str(numdeath) + " absorptions occurs during the simulation\n") 
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   ---------------------------FLUX(Z)----------------------------\n")  
    file.write(" --- Z ---         ---FLUX(Z)--- \n")
    z=0;
    for i in range(len(flux_local)):
        #file.write("   "+ str(round(z,4))+ " "*(6-len(str(round(z,2)))) + "          " + str(round(flux_local[i],6))+" "*(12-len(str(round(z,4)))) +str(round(test2[i],6)) +"\n");
        z = z + step
    file.close(); 
    for i in range(len(flux_local)):
        fluxPSP1.write("   "+ str(flux_local[i])+"\n");
    for i in range(len(flux_local)):
        coordPSP1.write("   "+ str(z)+"\n");
        z = z + step
    for i in range(len(flux_local)):
        variancePSP1.write("   "+ str(variancepic[i])+"\n");
    #for i in range(len(flux_local)):
        #deviationPSP1.write("   "+ str(variancepic[i])+"\n");
    deviationPSP1.close()                     
    variancePSP1.close()     
    fluxPSP1.close();
    coordPSP1.close();          
    return flux_local
             
def fsp(x,xi):                          # function xi - (1+z)exp(-z)
    return  xi - (1+x)*np.exp(-x)             

listC = [0.1,0.2,0.5,0.8,0.9,0.95,0.99,0.999]
for i in range(len(listC)):
    a = main(listC[i]);
    filename = str(listC[i])
    file = open(filename + '.txt', 'w') 
    file.write('max of the flux = ' +  str(max(a)) + '\n')
    file.write("average on bins max = "+  str((a[a.index(max(a))]+a[a.index(max(a))+1])/2))
    file.close()
            
print("---------------------------------Program-------Ends-----------------------------")