# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:42:57 2017

@author: TH Lab
"""
import random

def update_collision(z,w,c,z1,step,flux_local,iflux,iscat,ideath):
    index = int(z1/step)+1;                             #Count in the bins
    flux_local[index-1] = flux_local[index-1] + 1;      #Increment event for flux
    iflux[index-1] = iflux[index-1] +1;                 #Increment event for flux
    xi = random.random();
    if (xi < c):                                        #Check if it's a scattering with a probability of Sigmascat/sigmatot
        iscat = iscat + 1;                              #Increment of the number of scattering for the particle
        z=z1;                                           #Save the new position
        randnum=random.random()                         #Isotropic collision, direction random
        w=2*randnum-1;
    else:                                                #Check if it's an absorption with a probability Sigmaabs/sigmatot
        ideath = 1;                                      #Kill the particle
    return flux_local,iflux,iscat,z,w,ideath