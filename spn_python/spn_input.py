import numpy as np
from sys import exit

r"""
This class handles the user inputs:
- checks for errors and adapts inputs when applicable
- generates the array with spatial points for the solution and for source
- starts report file
"""

class spn_input(): #input material parameters

    def __init__(self, classical, number_of_points_in_solution,
            total_length_system, left_boundary_condition, right_boundary_condition,
            total_cross_section, scattering_ratio, source, moments_of_pdf,
            source_region, source_width, solver,filename):

        self.filename = filename

        # Check if input can be converted
        # Guarantees non-negative cross section
        try:
            self.xs = float(total_cross_section)
        except:
            print('\n Error in the total cross section \n')
            exit(0)
        if self.xs<0:
            self.xs = 0.

        # Check if nonclassical input of moments can be converted
        # Creates classical input of moments
        self.moments = moments_of_pdf
        if str(classical)!='y':
            for idx, val in enumerate(self.moments):
                try:
                    self.moments[idx] = float(val)
                except:
                    print('\n Error in the input of moments \n')
                    exit(0)
            self.classical = 0
        else:
            for i in range(6):
                self.moments[i] = np.math.factorial(i+1)/(self.xs**(i+1))
            self.classical = 1

        # Check if input can be converted
        # Guarantees minimum number of cells to be 10 (11 points)
        try:
            self.M = int(number_of_points_in_solution)
        except:
            print('\n Error in the number mesh cells \n')
            exit(0)
        if self.M < 10:
            self.M = 10

        # Check if input can be converted
        # If non-positive total length is entered, assumes 1
        try:
            self.length = float(total_length_system)
        except:
            print('\n Error in the the total length of the system\n')
            exit(0)
        if self.length<=0:
            self.length = 1.0

        # Check if input can be converted
        # Guarantees scattering ratio is 0 <= c <= 1
        try:
            self.c = float(scattering_ratio)
        except:
            print('\n Error in the scattering ratio\n')
            exit(0)
        if self.c<0:
            self.c = 0.
        elif self.c>1.0:
            self.c = 1.0

        # Check if input can be converted
        # Guarantees non-negative source
        try:
            self.source = float(source)
        except:
            print('\n Error in the source\n')
            exit(0)
        if self.source<0.:
            self.source=0.

        # Check if input can be converted
        # Guarantees non-negative left boundary condition
        try:
            self.bcl = float(left_boundary_condition)
        except:
            print('\n Error in the left boundary condition\n')
            exit(0)
        if self.bcl < 0.:
            self.bcl = 0.

        # Check if input can be converted
        # Guarantees non-negative right boundary condition
        try:
            self.bcr = float(right_boundary_condition)
        except:
            print('\n Error in the right boundary condition\n')
            exit(0)
        if self.bcr < 0.:
            self.bcr = 0.

        # Check if choice of solver is valid
        if solver not in ['bicg','bicgstab','cg','cgs','gmres','lgmres','minres','qmr']:
            print('\n Error in the choice of solver: do not know %r \n')
            exit(0)
        else:
            self.solver = solver

        # Check if input can be converted
        # Check that source_width is smaller than length and positive
        if str(source_region) is 'y':
            try:
                self.source_width = float(source_width)
            except:
                print('\n Error in the width of the source region\n')
                exit(0)
            if self.source_width >= self.length or self.source_width <= 0:
                print('\n Error in the width of the source region\n')
                exit(0)
            self.source_region=1
        else:
            self.source_region=0

    #    self.sca = self.c * self.xs # Scattering cross section
        self.abs = (1-self.c) * self.xs # Absorbing cross section
        self.x = np.linspace(0.,self.length,self.M+1) # Array of points where solution will be given
        self.source_array = np.full((self.M+1),self.source) # Full array of source

        # Adjust source array if source is only in center region
        if self.source_region==1:
            for idx,val in enumerate(self.x):
                if val <= self.length/2. - self.source_width/2 or val >= self.length/2. + self.source_width/2:
                    self.source_array[idx] = 0.

#        # Directions and weights for quadrature
#        if self.slab is 0: # If Rod Geometry
#            (self.mu,self.weight) = (np.array([1,-1]),np.array([1,1]))
#            self.weight = np.array([1,1])
#        else: # If slab geometry: Gauss-Legendre
#            (self.mu,self.weight) = np.polynomial.legendre.leggauss(self.N)
#            (self.mu,self.weight) = (self.mu[::-1],self.weight[::-1]) # Start at mu=1

        with open(self.filename+"n_report.txt","w") as file:
            file.write("*******************************************************\n")
            file.write("List of system and material parameters:\n\n")
            file.write("            Model: CLASSICAL\n\n" if self.classical==1
                        else "            Model: NON-CLASSICAL\n\n")
            if self.classical==1:
                file.write("         Total cross section: %r\n" %self.xs)
                file.write("     Absorbing cross section: %r\n" %self.abs)
            else:
                for idx,val in enumerate(self.moments):
                    file.write("                    Moment %r: " %(idx+1) + "%r\n" %val)
            file.write("            Scattering ratio: %r\n" %self.c)
            file.write("                      Source: %r\n" %self.source)
            file.write("               Source region: %r\n\n" %self.source_width if self.source_region==1 else "               (No source region)\n\n")
            file.write("  Total length of the system: %r\n" %self.length)
            file.write("     Left boundary condition: %r\n" %self.bcl)
            file.write("    Right boundary condition: %r\n\n" %self.bcr)
            file.write("     Number of cells in mesh: %r\n" %self.M)
            file.write("            Iterative solver: %r\n" %self.solver)
            file.write("                Name of file: %r\n" %self.filename)
            file.write("*******************************************************\n\n")
            file.write("*******************************************************\n")
