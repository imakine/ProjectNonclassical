# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 12:38:42 2017

@author: TH Lab
"""
import random
import moment as m
import check_out as check
import update_collision as u
import scipy.optimize
import numpy as np

def lifeOfParticle(sp,iesc,ideath,sigmatot,freepath,s,z,w,c,step,flux_local,iflux,iscat,thickness):
    while((iesc+ideath)==0):
        if (sigmatot !=0.):
            randnum=random.random();
            if (sp==1):
                freepath = scipy.optimize.brentq(fsp,0.0,100.0,args=(randnum))/((3**0.5)*sigmatot) #Sampling of the distance travelled
            elif (sp==2):
                if randnum < 4/9:
                    freepath = 0;
                else:
                    xi = 9*(1-randnum)/5   # adjust random number xi to xi_2 such that: 0 < xi_2 < 1
                    freepath = scipy.optimize.brentq(fsp,0.0,100.0,args=(xi))*(3/5)**0.5/sigmatot
            elif (sp==3):
                try:
                    freepath = scipy.optimize.brentq(f,0,100,args=((randnum)))/sigmatot
                except:
                    print('module scipy non efficient')
                    freepath = secante(0,100,1e-6,randnum,sigmatot)
            s = m.moment(s,freepath);
        z1=z+w*freepath;                                            #Update of the position
        iesc =check.check_out(z1,thickness,sigmatot);
        if (not iesc):
            flux_local,iflux,iscat,z,w,ideath = u.update_collision(z,w,c,z1,step,flux_local,iflux,iscat,ideath)
    return s,iesc,ideath,iscat,flux_local,iflux 

def fsp(x,xi):                          # function xi - (1+z)exp(-z)
    return  xi - (1+x)*np.exp(-x)

def f(x,xi):
    a = 5.642025;
    b = 0.469086;
    l = (5+2*(10/3)**(0.5))**(0.5);
    m = (5-2*(10/3)**(0.5))**(0.5);
    y = a/(l**2);
    z = b/(m**2);
    return (xi - y*(1-fsp(l*x)) - z*(1-fsp(m*x)))

def secante(a,b,prec,xi,sigmatot):
	while f(a*sigmatot,xi) > prec:
		a = a-f(a*sigmatot,xi)*(b-a)/(f(b*sigmatot,xi)-f(a*sigmatot,xi))
	return a 