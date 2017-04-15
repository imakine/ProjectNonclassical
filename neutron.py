import random
import math

def randomGen():
    xi = random.randrange(1,1001)
    xi = xi / 1000
    return xi

def main():

    #Initialisation du generateur de nombres aleatroires
    randnum=randomGen();
    print("1st Random Number : %f\n",randnum);

    #Initialisation des parametree de la simulation

    #Blindage
    thickness = float(input("Thickness of the blindage : "));

    #Cross-sections
    sigmaabs = float(input("Absorption cross-section : "));
    sigmascat = float(input("Scattering cross-section : "));
    sigmatot=sigmaabs+sigmascat;
    lamb = 1/sigmatot; #mean-free path

    #Number of path

    n = int(input("Nombre de trajectoires : "));


    # Parametres particule incidente

    x0=0;
    y0=0;
    z0=0;
    theta = float(input("Angle d'incidence : "));
    u0=math.sin(math.pi/180*theta);
    v0=0;
    w0=math.cos(math.pi/180*theta);

    #Initialisation nombres particules perdues, absorbees ou emises */
    global numlost;
    global numabs;
    global numesc;
    global numdeath;
    numlost=0;
    numabs=0;
    numesc=0;
    numdeath =0;
    for i in range(1,n+1):
        x=x0;
        y=y0;
        z=z0;
        u=u0;
        v=v0;
        w=w0;
        ilost=0;
        iesc=0;
        ideath=0;
        weight= 1;
        while((ilost+iesc+ideath)==0): 
            randnum=randomGen();
            freepath=-lamb*(math.log(randnum));
            x1=x+u*freepath;
            y1=y+v*freepath;
            z1=z+w*freepath;
            if(z1 >= thickness):
                iesc=weight;
            else: 
                if(z1 <= 0.):
                    ilost=weight;
                else:
                    weight=weight*sigmascat/sigmatot;
                    x=x1;
                    y=y1;
                    z=z1;
                    randnum=randomGen();
                    phi=2*math.pi*randnum;
                    randnum=randomGen();
                    w=2*randnum-1;
                    thet=math.acos(w);
                    u=math.sin(thet)*math.cos(phi);
                    v=math.sin(thet)*math.sin(phi);
                    if(weight<1e-6):
                        randnum =randomGen();
                        if(randnum <1/2):
                            ideath = 1;
                        else:
                            weight = 2*weight;
        numlost=numlost+ilost;
        #numabs+=iabs;
        numesc=numesc+iesc;
        numabs=numdeath+ideath;
    print("n =",n);
    print("numlost =",numlost);
    print("numabs =",numabs);
    print("numesc=",numesc);
