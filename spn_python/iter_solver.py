import scipy.sparse.linalg  as spla

class iter_solver():

    def __init__(self,matrix,vector,solver,info,filename):

        limit=1000
        tolerance = 1e-6
        if solver == 'bicg': #BIConjugate Gradient iteration
            (self.sol,self.info) = spla.bicg(matrix,vector,tol=tolerance,maxiter=limit)

        elif solver == 'bicgstab': #BIConjugate Gradient STABilized iteration
            (self.sol,self.info) = spla.bicgstab(matrix,vector,tol=tolerance,maxiter=limit)

        elif solver == 'cg': #Conjugate Gradient iteration
            (self.sol,self.info) = spla.cg(matrix,vector,tol=tolerance,maxiter=limit)

        elif solver == 'cgs': #Conjugate Gradient Squared iteration
            (self.sol,self.info) = spla.cgs(matrix,vector,tol=tolerance,maxiter=limit)

        elif solver == 'gmres': #Generalized Minimal RESidual iteration
            (self.sol,self.info) = spla.gmres(matrix,vector,tol=tolerance,maxiter=limit)

        elif solver == 'lgmres': #Loose Generalized Minimal RESidual iteration
            (self.sol,self.info) = spla.lgmres(matrix,vector,tol=tolerance,maxiter=limit)

        elif solver == 'minres': #MINimum RESidual iteration
            (self.sol,self.info) = spla.minres(matrix,vector,tol=tolerance,maxiter=limit)

        elif solver == 'qmr': #Quasi-Minimal Residual iteration
            (self.sol,self.info) = spla.qmr(matrix,vector,tol=tolerance,maxiter=limit)

        with open(filename+"n_report.txt","a") as file:
            file.write("Solver %r was SUCCESSFUL in " %solver + info + "\n" if self.info==0
            else "ATTENTION\n")
            if self.info !=0:
                file.write("Convergence to tolerance not achieved after %r iterations in "
                    %self.info + "solver %r for " %solver + info + "\n" if self.info>0
                    else "Solver %r presents an illegal input or breakdown in " %solver + info + "\n")

#Provides convergence information:
#0 : successful exit >0 : convergence to tolerance not achieved, number of iterations <0 : illegal input or breakdown
