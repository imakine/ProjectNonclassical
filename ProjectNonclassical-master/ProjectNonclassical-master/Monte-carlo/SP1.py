import random
import math
import time
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize 


def randomGen():
    xi = random.random()
    xi = float(xi);
    return xi

def main():

    print("1D code is running with theta AND phi");

    #Initialisation des parametree de la simulation

    #Blindage
    thickness = int(input("Thickness of the blindage (integer form) : "));

    #Cross-sections
    sigmaabs = float(input("Absorption cross-section [thickness^-1] : "));
    sigmascat = float(input("Scattering cross-section [thickness^-1] : "));
    sigmatot=sigmaabs+sigmascat;
    global lamb;
    if (sigmatot !=0):
        lamb = 2/((3)**(0.5)*sigmatot); #mean-free path
        print("lamb",lamb)
    #Source (discrete)

    Q = float(input("Difusion source (discrete) [neutrons.thickness^-1.s^-1] = "));
    
    #Number of path

    n = int(input("Number of path : "));

    #step

    step = float(input("Number of step for flux: (0.01, 0.1, 1) = "))

    #Output file

    filename = str(input("Name of the output file : "));
    file = open(filename + ".txt", "w") 
    fluxSP1 = open("fluxSP1.txt","w")
    varianceSP1 = open("varianceSP1.txt","w");
    deviationSP1 = open("deviationSP1.txt","w");
    coordSP1 = open("coordSP1.txt","w")
    
    # Parametres particule incidente

    initialtime = time.time();
   

    #Initialisation nombres particules perdues, absorbees ou emises */
    global numlost;
    global numabs;
    global numesc;
    global numdeath;
    global numscattot;
    global flux_col;
    global flux_abs;
    global flux_track;
    global freepath;
    global flux_local
    global s
    global bybin;
    global variancepic;
    global deviationpic;
    test = test2 = fluxpic= variancepic = [0]*(int(thickness/step));
    s =[0]*6;
    flux_local = [0]*(int(thickness/step)); # collision between 0 and 1 is in first place, between 1 and 2 is second ...
    numlost=0;
    numabs=0;
    numesc=0;
    numdeath =0;
    numscattot = 0;
    flux_col = 0;
    flux_abs =0;
    flux_track =0;
    freepath =0;
    for i in range(1,n+1):
         #x0=0;
        #y0=0;
        randnum = randomGen();
        z0=thickness*randnum; #thickness
        randnum = randomGen();
        theta = 180*randnum;
        #u0=math.sin(math.pi/180*theta);
        #v0=0;
        w0=math.cos(math.pi/180*theta);
        #x=x0;
        #y=y0;
        z=z0;
        #u=u0;
        #v=v0;
        w=w0;
        ilost=0;
        iesc=0;
        ideath=0;
        iscat=0;
        iflux_col = 0;
        iflux_abs =0 ;
        iflux_track =0;
        ivariancepic = [0]*(int(thickness/step));
        iflux = [0]*(int(thickness/step));
        #weight= 1;
        while((ilost+iesc+ideath)==0):
            randnum=randomGen();                
            if (sigmatot !=0):
                freepath = scipy.optimize.brentq(fsp,0.0,100.0,args=(randnum))/((3**0.5)*sigmatot)
                #freepath= secante(0,1,1e-7,randnum,sigmatot);
                s[0] = s[0] + freepath;
                s[1] = s[1] + freepath**2;
                s[2] = s[2] + freepath**3;
                s[3] = s[3] + freepath**4;
                s[4] = s[4] + freepath**5;
                s[5] = s[5] + freepath**6;
            #x1=x+u*freepath;
            #y1=y+v*freepath;
            z1=z+w*freepath;
            if(z1 >= thickness or sigmatot==0):
                iesc=1;
            elif(z1 < 0.):
                ilost=1; 
            else:
                xi = randomGen();
                if (xi < sigmascat/sigmatot): #that's a scattering
                    index = int(z1/step)+1;
                    flux_local[index-1] = flux_local[index-1] + 1;  
                    iflux[index-1] = iflux[index-1] +1;
                    ivariancepic[index-1] = ivariancepic[index-1] + 1;          
                    iscat = iscat + 1;
                    iflux_col = iflux_col + 1
                    iflux_track =iflux_track + freepath
                    #x=x1;
                    #y=y1;
                    z=z1;
                    #randnum=randomGen();
                    #phi=2*math.pi*randnum;
                    randnum=randomGen();
                    w=2*randnum-1;
                    #thet=math.acos(w);
                    #u=math.sqrt(thet)*math.cos(phi);
                    #v=math.sin(thet)*math.sin(phi);
                if (sigmaabs !=0 and xi > sigmascat/sigmatot): #that's an absorption
                    index = int(z1/step) + 1;
                    flux_local[index-1] = flux_local[index-1] + 1; 
                    iflux[index-1] = iflux[index-1] +1;
                    ivariancepic[index-1] = ivariancepic[index-1] + 1;          
                    ideath = 1; 
                    iflux_abs = iflux_abs + 1; 
                else:
                    ideath = 0; 
        for i in range(len(fluxpic)):            
            fluxpic[i] = fluxpic[i] + iflux[i];            
            variancepic[i] = variancepic[i] + ivariancepic[i]**2;                      
        flux_track = flux_track + iflux_track;                    
        flux_abs = flux_abs + iflux_abs;                    
        flux_col = flux_col + iflux_col;                    
        numlost=numlost+ilost;
        numesc=numesc+iesc;
        numdeath=numdeath+ideath;
        numscattot = numscattot + iscat;
    flux_track = flux_track/n;    
    flux_col = flux_col/n;
    flux_abs = flux_abs/n;
    for i in range(len(fluxpic)):
        variancepic[i] = variancepic[i]/(n-1)*((Q*thickness)/step/sigmatot)**2;
        fluxpic[i] = fluxpic[i]/n*(Q*thickness)/step/sigmatot;                                                  
        test[i] = variancepic[i]-(fluxpic[i])**2;
        test2[i] = test[i]**(0.5);
    flux_local[:] = [(x/n)*(Q*thickness)/step/sigmatot for x in flux_local];
    s[:] = [x / (numscattot+numdeath+numesc+numlost) for x in s];        
    finaltime = time.time();
    lon = thickness//2;
    source = sigmaabs*(flux_local[lon-2]+flux_local[lon-1] + flux_local[lon] +flux_local[lon+1]+flux_local[lon+2])/5;
    print("n =",n);
    print("numscattot = ", numscattot);
    print("numlost =",numlost);
    print("numabs =",numdeath);
    print("numesc=",numesc);
    print("flux_col = " ,flux_col);
    print("flux_abs = " ,flux_abs);
    print("flux_track = " ,flux_track);
    print("flux(z) = ",flux_local);
    print("s = ",s);   
    print("source = ",source);
    print("Time elapsed during the running of the code : ",finaltime - initialtime, "seconds");
    print(" ");
    print("Variance of flux = ",test)  
    print(" ");
    print("Standart deviation of flux = ", test2)
    print(" ");
    print("Maximum for the standart deviation = "+ str(max(test2)) + " at the position " + str(test2.index(max(test2))) )
    plt.plot(flux_local)
    plt.ylabel('Flux(z)')
    plt.show()
    file.write("\n")
    file.write("\n")
    file.write("\n")
    file.write("   ------------------------------------------------------------\n")
    file.write("   |                                                          |\n")
    file.write("   |                 MONTE-CARLO CODE FOR                     |\n")
    file.write("   |              SP1 EQUATION (1D- Continuous source)        |\n")
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
        file.write("   "+ str(round(z,4))+ " "*(6-len(str(round(z,2)))) + "          " + str(round(flux_local[i],6))+" "*(12-len(str(round(z,4)))) +str(round(test2[i],6)) +"\n");
        z = z + step
    file.close(); 
    for i in range(len(flux_local)):
        fluxSP1.write("   "+ str(flux_local[i])+"\n");
    for i in range(len(flux_local)):
        coordSP1.write("   "+ str(z)+"\n");
        z = z + step
    for i in range(len(flux_local)):
        varianceSP1.write("   "+ str(variancepic[i])+"\n");
    for i in range(len(flux_local)):
        deviationSP1.write("   "+ str(variancepic[i])+"\n");
    deviationSP1.close()                     
    varianceSP1.close()    
    fluxSP1.close();
    coordSP1.close();                 
    return flux_local

def col(liste):
    for i in range(len(liste)):
        split_num = str(liste[i]).split('.');
        int_part = str(split_num[0]);
        decimal_part = str(split_num[1]);
        print(int_part + "," + decimal_part);
             
def fsp(x,xi): # function xi - (1+z)exp(-z)
    return  xi - (1+x)*np.exp(-x)             
             
             
def f(s,xi,sigmatot):
    return xi - (1-(1+(3)**(0.5)*sigmatot*s)*(math.e)**(-(3)**(0.5)*sigmatot*s))

def secante(a,b,prec,xi,sigmatot):
    while f(a,xi,sigmatot) > prec:
        
        c=a;
        a = a-f(a,xi,sigmatot)*(b-a)/(f(b,xi,sigmatot)-f(a,xi,sigmatot))
        b = c;  
    return a

        
            
print("---------------------------------Program-------Ends-----------------------------")
