import numpy as np
import scipy.linalg
from iter_solver import *


class spn_1():

    def __init__(self,inp):

        # Mesh interval
        self.h = inp.length/inp.M
        # Modified Sigma_t
        self.Et = 2*inp.moments[0]/inp.moments[1]
        # Modified Sigma_a
        self.Ea = (1-inp.c)/inp.moments[0]
        # Diffusion coefficient
        self.D = 1/(3*self.Et)

        # Create matrix A
        self.diag = np.full((inp.M+1),(self.Ea+2*self.D/(self.h**2))) # Diagonal terms
        self.offdiag = np.full(inp.M,(-self.D/(self.h**2))) # Upper and lower diagonal terms
        self.A = np.diag(self.offdiag,-1) + np.diag(self.diag) + np.diag(self.offdiag,1) # create matrix

        # Adjust first and last rows of matrix A for boundary conditions
        self.A[0,0] = 1./2. + 3*self.D/(2*self.h)
        self.A[0,1] = -4*self.D/(2*self.h)
        self.A[0,2] = self.D/(2*self.h)
        self.A[inp.M,inp.M] = self.A[0,0]
        self.A[inp.M,(inp.M-1)] = self.A[0,1]
        self.A[inp.M,(inp.M-2)] = self.A[0,2]

        # Solve system
        self.iter = iter_solver(self.A,inp.source_array,inp.solver,'solving SP1 system',inp.filename)
        self.scalar = self.iter.sol
        
