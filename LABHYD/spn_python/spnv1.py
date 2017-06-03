from spn_main import *
import numpy as np
#***********************************************************
#***********************************************************
# User inputs
filename = 'testout' # Name of files generated
solver = 'lgmres' # Iterative solver to be used in 'iter_solver'

classical = 'y'    # Classical SPN (y/n)? (Anything besides 'y' defines nonclassical)

M = 12500   # Number of cells in mesh (mesh intervals)

length = 100    # Total length of system
bcl = 0. #'r'   # Left incident boundary condition
bcr = 0. #'r'   # Right incident boundary condition

xs = 1.0    # Total cross section
c = 0.8  # Scattering ratio
Q = 1    # Homogeneous source

# Moments of the pdf (input in case of nonclassical)
moments = np.zeros(6)
moments[0] = 1.  # First moment of the pdf
moments[1] = 1.  # Second moment of the pdf
moments[2] = 1.  # Third ...
moments[3] = 1.  #
moments[4] = 1.  #
moments[5] = 1.  #

source_region = 'y' # Is the source a region around the center (y/n)? (Anything besides 'y' will define the source everywhere)
# If source_region = 'y', enter source_width below
source_width = 1 # TOTAL width of the source around the center (eg. 1.0 means source is in (center-0.5, center+0.5)

#***********************************************************
#Starts time counter
code_time = time.time()
#***********************************************************

# Runs code
spn_main(classical,M,length,bcl,bcr,xs,c,Q,moments,source_region,source_width,solver,code_time,filename+'_sp') # Runs code
