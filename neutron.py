import random
import math

def randomGen():
    xi = random.randrange(1,1001)
    xi = xi / 1000
    xi = float(xi);
    return xi

def main():

    print("1D code is running with theta AND phi");

    #Initialisation des parametree de la simulation

    #Blindage
    thickness = float(input("Thickness of the blindage [cm] : "));

    #Cross-sections
    sigmaabs = float(input("Absorption cross-section [cm^-1] : "));
    sigmascat = float(input("Scattering cross-section [cm^-1] : "));
    sigmatot=sigmaabs+sigmascat;
    lamb = 1/sigmatot; #mean-free path

    #Number of path

    n = int(input("Number of path : "));


    # Parametres particule incidente

    x0=0;
    y0=0;
    z0=0;
    theta = float(input("Initial angle [°] : "));
    u0=math.sin(math.pi/180*theta);
    v0=0;
    w0=math.cos(math.pi/180*theta);

    #Initialisation nombres particules perdues, absorbees ou emises */
    global numlost;
    global numabs;
    global numesc;
    global numdeath;
    global numscattot;
    global flux_col;
    global flux_abs;
    global flux_track;
    numlost=0;
    numabs=0;
    numesc=0;
    numdeath =0;
    numscattot = 0;
    flux_col = 0;
    flux_abs =0;
    flux_track =0;
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
        iscat=0;
        iflux_col = 0;
        iflux_abs =0 ;
        iflux_track =0;
        weight= 1;
        while((ilost+iesc+ideath)==0): 
            randnum=randomGen();
            freepath=-lamb*(math.log(randnum));
            x1=x+u*freepath;
            y1=y+v*freepath;
            z1=z+w*freepath;
            if(z1 >= thickness):
                iesc=1; #weight ?
            if(z1 <= 0.):
                ilost=1; #weight ?
            else:
                iscat = iscat + 1;
                iflux_col = iflux_col + weight*(1/sigmatot)
                iflux_track =iflux_track + weight*freepath
                weight=weight*sigmascat/sigmatot;
                xi = randomGen();
                if (xi < sigmascat/sigmatot): #that's a scattering
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
                    #if(weight<1e-6):
                     #   randnum =randomGen();
                     #   if(randnum <1/2):
                      #      ideath = 1;
                      #      iflux_abs = iflux_abs + weight*(1/sigmaabs);
                       # else:
                       #     weight = 2*weight;
                if (sigmaabs !=0 and xi > sigmascat/sigmatot): #that's an absorption
                    ideath = 1;
                    iflux_abs = iflux_abs + weight*(1/sigmaabs);
                else:
                    ideath = 0;
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
    print("n =",n);
    print("numscattot = ", numscattot);
    print("numlost =",numlost);
    print("numabs =",numdeath);
    print("numesc=",numesc);
    print("flux_col = " ,flux_col);
    print("flux_abs = " ,flux_abs);
    print("flux_track = " ,flux_track);
    print("---------------------------------Program-------Ends-----------------------------")
