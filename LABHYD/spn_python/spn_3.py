import numpy as np
import scipy.linalg
from iter_solver import *


class spn_3():

    def __init__(self,inp):

        # Mesh interval
        self.h = inp.length/inp.M
        # Lambda_1,Lambda_2
        self.lambda1 = 3*inp.moments[3]/(10*inp.moments[1]**2)-inp.moments[2]/(3*inp.moments[0]*inp.moments[1])
        self.lambda2 = (9*inp.moments[4]/5-27*inp.moments[0]*inp.moments[5]/(21*inp.moments[1])+
            3*inp.moments[2]*inp.moments[3]/inp.moments[1]-10*inp.moments[2]**2/(3*inp.moments[0]
            ))/(10*inp.moments[1]*inp.moments[2]-9*inp.moments[0]*inp.moments[3])
        # Beta_1,Beta_2
        self.beta1 = inp.moments[2]/(3*inp.moments[0]*inp.moments[1])-1
        self.beta2 = (10*inp.moments[2]**2/(3*inp.moments[0])-9*inp.moments[4]/5
            )/(10*inp.moments[1]*inp.moments[2]-9*inp.moments[0]*inp.moments[3])-1

        # Modified Sigma_t
        self.Et = 2*inp.moments[0]/inp.moments[1]
        # Modified Sigma_a
        self.Ea = ((1-inp.c)/inp.moments[0])/(1+self.beta1*(1-inp.c))
        # Modified Source
        self.modsource = inp.source_array/(1+self.beta1*(1-inp.c))
        # Modified Sigma_2
        self.E2 = 4*(1+self.beta1*(1-inp.c))*(1-self.beta2*(1-inp.c))/(5*self.lambda1*inp.moments[0])
        # Modified Sigma_3
        self.E3 = 27*self.lambda1*self.Et/(28*self.lambda2*(1+self.beta1*(1-inp.c))-self.lambda1)
        # Diffusion coefficient
        self.D = 1/(3*self.Et)


        # Create block Matrix A_11
        self.diag11 = np.full((inp.M+1),(self.Ea+2*self.D/(self.h**2))) # Diagonal terms
        self.offdiag11 = np.full(inp.M,(-self.D/(self.h**2))) # Upper and lower diagonal terms
        self.A11 = np.diag(self.offdiag11,-1) + np.diag(self.diag11) + np.diag(self.offdiag11,1)
        # Create block Matrix A_12
        self.diag12 = np.full((inp.M+1),(4*self.D/(self.h**2)))
        self.offdiag12 = 2*self.offdiag11
        self.A12 = np.diag(self.offdiag12,-1) + np.diag(self.diag12) + np.diag(self.offdiag12,1)
        # Create block Matrix A_21
        self.diag21 = np.full((inp.M+1),(4*self.D/(5*self.h**2)))
        self.offdiag21 = np.full(inp.M,(-2*self.D/(5*self.h**2)))
        self.A21 = np.diag(self.offdiag21,-1) + np.diag(self.diag21) + np.diag(self.offdiag21,1)
        # Create block Matrix A_22
        self.diag22 = np.full((inp.M+1),(self.E2+(2*self.D/(self.h**2))*(4/5+27*self.Et/(35*self.E3))))
        self.offdiag22 = np.full(inp.M,(-(4/5+27*self.Et/(35*self.E3))*self.D/(self.h**2)))
        self.A22 = np.diag(self.offdiag22,-1) + np.diag(self.diag22) + np.diag(self.offdiag22,1)

        self.A = np.bmat([[self.A11, self.A12], [self.A21, self.A22]]) # Create matrix A

        # Adjust rows of matrix A for left boundary condition
        self.A[0,0] = 1./2. + 3*self.D/(2*self.h)
        self.A[0,1] = -4*self.D/(2*self.h)
        self.A[0,2] = self.D/(2*self.h)
        self.A[0,(inp.M+1)] = 5./8 + 6*self.D/(2*self.h)
        self.A[0,(inp.M+2)] = -8*self.D/(2*self.h)
        self.A[0,(inp.M+3)] = 2*self.D/(2*self.h)
        self.A[(inp.M+1),0] = -1./8
        self.A[(inp.M+1),1] = 0.
        self.A[(inp.M+1),2] = 0.
        self.A[(inp.M+1),(inp.M+1)] = 5./8 + (3/(7*self.E3))*3/(2*self.h)
        self.A[(inp.M+1),(inp.M+2)] = -(3/(7*self.E3))*4/(2*self.h)
        self.A[(inp.M+1),(inp.M+3)] = (3/(7*self.E3))*1/(2*self.h)

        # Adjust rows of matrix A for right boundary condition
        self.A[inp.M,inp.M] = self.A[0,0]
        self.A[inp.M,(inp.M-1)] = self.A[0,1]
        self.A[inp.M,(inp.M-2)] = self.A[0,2]
        self.A[inp.M,(2*inp.M+1)] = self.A[0,(inp.M+1)]
        self.A[inp.M,(2*inp.M)] = self.A[0,(inp.M+2)]
        self.A[inp.M,(2*inp.M-1)] = self.A[0,(inp.M+3)]
        self.A[(2*inp.M+1),inp.M] = self.A[(inp.M+1),0]
        self.A[(2*inp.M+1),(inp.M-1)] = self.A[(inp.M+1),1]
        self.A[(2*inp.M+1),(inp.M-2)] = self.A[(inp.M+1),2]
        self.A[(2*inp.M+1),(2*inp.M+1)] = self.A[(inp.M+1),(inp.M+1)]
        self.A[(2*inp.M+1),(2*inp.M)] = self.A[(inp.M+1),(inp.M+2)]
        self.A[(2*inp.M+1),(2*inp.M-1)] = self.A[(inp.M+1),(inp.M+3)]

        # Create vector B
        self.B = np.zeros((2*inp.M+2))
        self.B[0:(inp.M+1)] = self.modsource
        # Solve system
        self.iter = iter_solver(self.A,self.B,inp.solver,'solving SP3 system',inp.filename)
        self.scalar = self.iter.sol[0:(inp.M+1)]
        