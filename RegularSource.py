import random
import math
import time
 

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
        lamb = 1/sigmatot; #mean-free path

    #Number of path

    n = int(input("Number of path : "));


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
    s =[0]*6;
    flux_local = [0]*(thickness+1); # collision between 0 and 1 is in first place, between 1 and 2 is second ...
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
        ilam = 0;
        #weight= 1;
        while((ilost+iesc+ideath)==0):
            randnum=randomGen();
            if (sigmatot !=0):
                freepath=-lamb*(math.log(randnum));
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
                    index = int(z1) + 1;
                    flux_local[index] = flux_local[index] + 1;
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
                    ideath = 1; 
                    iflux_abs = iflux_abs + 1; 
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
    flux_local[:] = [x/n for x in flux_local]; # x/n
    if (numscattot !=0):
        s[:] = [x / numscattot for x in s];
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
    return flux_local;
def col(liste):
    for i in range(len(liste)):
        split_num = str(liste[i]).split('.');
        int_part = str(split_num[0]);
        decimal_part = str(split_num[1]);
        print(int_part + "," + decimal_part);
        
            
    print("---------------------------------Program-------Ends-----------------------------")
