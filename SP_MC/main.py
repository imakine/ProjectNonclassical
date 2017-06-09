# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:44:29 2017

@author: TH Lab
"""
from Initialisation import initialisation
from initializeOneParticle import initializeOneParticle
from lifeOfParticle import lifeOfParticle
from addParticleInSample import addParticleInSample
from averageSample import averageSample
import time
from output import output
from theoretical import theoretical

def main(c,thickness):
    '''
    This function gives the flux for SP1 or SP2 or SP3 with Monte-Carlo Simulation
    c is the scattering ratio (ratio between the scattering cross section and the total cross-section)
    '''
    
    #Initialisation of the parameters of the simulation
    
    #thickness = 80;                     #Thickness of the mediq
    sigmatot = 1;                       #Total cross section
    sigmaabs = (1-c)*sigmatot;          #Absorption cross-section of the media
    sigmascat = c*sigmatot;             #Scattering cross-section of the media
    region = True                       # If True the source is regional source on a distance sourceregion about the half of the media, if False, the source is uniform
    sourceregion = 1.0                  #thickness of center source region    
    Q = 1                               #Source of neutrons
    n = 1000000;                       #Number of path
    step = 0.1;                         #Step for bins Number of bins will be Length/step
    sp = 2;                             #SP =1,2,3 
    #classical = True;                  #Classical T/f
    

    # Initialisation in the Memory.
    
    initialtime,variance,std,ERRORpic,fluxpic,variancepic,s,flux_local,numdeath,numscattot,numesc,freepath = initialisation(thickness,step)

    # Loop on each particle of the sample
    
    for i in range(1,n+1):
        
        # Interface with user
        
        if (i%1000 == 0):
            print(str(i/n*100) + "% of the running code done")  
            
        # Initialisation of the position, direction, etc. of the particle i
        
        z,w,iesc,ideath,iscat,iflux = initializeOneParticle(thickness,sourceregion,step,region)
        
        # Following of the life of the particle i
        
        s,iesc,ideath,iscat,flux_local,iflux=lifeOfParticle(sp,iesc,ideath,sigmatot,freepath,s,z,w,c,step,flux_local,iflux,iscat,thickness)
        
        #Computation of interest elements for the sample
        
        fluxpic,variancepic,numesc,numdeath,numscattot = addParticleInSample(flux_local,fluxpic,iflux,variancepic,numesc,iesc,numdeath,numscattot,iscat,ideath)
    fluxpic,variancepic,std,ERRORpic,flux_local,s = averageSample(fluxpic,variancepic,std,ERRORpic,n,Q,step,sigmatot,flux_local,s,variance,numscattot,numdeath,numesc)
    finaltime = time.time(); 
    
    # Computation of the different theoretical moments to obtain
    
    s1,s2,s3,s4,s5,s6 = theoretical(sp,sigmatot)
    
    # Generate all the output files
    
    output(s1,s2,s3,s4,s5,s6,s,finaltime,initialtime,flux_local,thickness,sigmaabs,sigmascat,sigmatot,Q,n,step,numscattot,numesc,numdeath,ERRORpic,std,variancepic)
    return flux_local

c = [0.95,0.99,0.999]
thickness = [100,200,600]
for i in range(2):
    a = main(c[i],thickness[i])
    filename = str("TEST")
    file = open(filename + str(c[i]) + '.txt', 'w')
    file.write('max of the flux = ' +  str(max(a)) + '\n')
    file.write("average on bins max = "+  str((a[a.index(max(a))]+a[a.index(max(a))+1])/2))
    file.close()