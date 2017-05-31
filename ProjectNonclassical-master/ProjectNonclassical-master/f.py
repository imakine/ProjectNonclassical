import random
import math

def classical():
    xi = random.randrange(0,1001)
    xi = xi / 1000
    return xi

def diffusion():
    sigmatot = 1;
    xi = classical();
    s = - math.log(1-xi);
    

