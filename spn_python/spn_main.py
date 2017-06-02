from spn_input import *
from spn_1 import *
from spn_2 import *
from spn_3 import *
from spn_output import *
import numpy as np
import time

r"""
This is the main file to solve SPn
"""

class spn_main():

    def __init__(self, classical, number_of_points_in_solution,
            total_length_system, left_boundary_condition, right_boundary_condition,
            total_cross_section, scattering_ratio, source, moments_of_pdf,
            source_region, source_width, solver,code_time,filename):

        # Class: handles user inputs, creates output "report"
        main_input = spn_input(classical, number_of_points_in_solution,
                total_length_system, left_boundary_condition, right_boundary_condition,
                total_cross_section, scattering_ratio, source, moments_of_pdf,
                source_region, source_width, solver,filename)

        # Class: solve SP1 system, calculates fluxes
        main_sp1 = spn_1(main_input)

        # Class: solve SP2 system, calculates fluxes
        main_sp2 = spn_2(main_input)

        # Class: solve SP3 system, calculates fluxes
        main_sp3 = spn_3(main_input)

        #Class: handles output of the code
        main_output = spn_output(main_input,main_sp1,main_sp2,main_sp3,code_time)
