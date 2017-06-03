import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time

class spn_output():

    def __init__(self, inp, sp1,sp2,code_time):

        #Ends time counter
        code_time = time.time() - code_time

        # Finishes "report" file with outputs and locations of other files
        with open(inp.filename+"n_report.txt","a") as file:
            file.write("*******************************************************\n\n")
            file.write("*******************************************************\n")
            file.write("List of outputs:\n\n")
            if inp.M%2==0:
                file.write("                       Center of the system: %r\n" %inp.x[int(inp.M/2)])
                file.write("SP1 scalar flux at the center of the system: %r\n" %sp1.scalar[int(inp.M/2)])
                file.write("SP2 scalar flux at the center of the system: %r\n" %sp2.scalar[int(inp.M/2)])
                #file.write("SP3 scalar flux at the center of the system: %r\n\n" %sp3.scalar[int(inp.M/2)])
            file.write("  Dimensions of scalar flux saved in " + inp.filename + "n_dim.csv\n")
            file.write("       Array of coordinates saved in " + inp.filename + "n_coord.csv\n\n")
            file.write("            SP1 scalar flux saved in " + inp.filename + "1_sclr.csv\n")
            file.write("    Plot of SP1 scalar flux saved in " + inp.filename + "1_plot.png\n\n")
            file.write("            SP2 scalar flux saved in " + inp.filename + "2_sclr.csv\n")
            file.write("    Plot of SP2 scalar flux saved in " + inp.filename + "2_plot.png\n\n")
            file.write("            SP3 scalar flux saved in " + inp.filename + "3_sclr.csv\n")
            file.write("    Plot of SP3 scalar flux saved in " + inp.filename + "3_plot.png\n\n")
            file.write("Plot of all SPn scalar flux saved in " + inp.filename + "n_plot.png\n")
            file.write("*******************************************************\n\n")
            file.write("*******************************************************\n")
            file.write("            Time elapsed: %r\n" %code_time)
            file.write("*******************************************************")

        # Creates csv file for dimensions (number of points)
        dimensions = np.array([inp.M+1])
        df = pd.DataFrame(dimensions)
        df.to_csv(inp.filename + 'n_dim.csv',header=None,index=False)

        # Creates csv file for Spatial coordinates
        df = pd.DataFrame(inp.x)
        df.to_csv(inp.filename + 'n_coord.csv',header=None,index=False)

        # Creates csv file for scalar flux
        df = pd.DataFrame(sp1.scalar)
        print('max sp1 = ' ,max(sp1.scalar))
        df.to_csv(inp.filename + '1_sclr.csv',header=None,index=False)
        df = pd.DataFrame(sp2.scalar)
        print('max sp2 = ' ,max(sp2.scalar))
        df.to_csv(inp.filename + '2_sclr.csv',header=None,index=False)
        #df = pd.DataFrame(sp3.scalar)
        #print('max sp3 = ' ,max(sp3.scalar))
        #df.to_csv(inp.filename + '3_sclr.csv',header=None,index=False)

        # Creates plot of scalar flux
        plt.plot(inp.x, sp1.scalar, 'b')
        plt.plot(inp.x, sp2.scalar, 'r')
        #plt.plot(inp.x, sp3.scalar, 'g')
        plt.savefig(inp.filename + 'n_plot.png')
        plt.clf()
        plt.plot(inp.x, sp1.scalar, 'b')
        plt.savefig(inp.filename + '1_plot.png')
        plt.clf()
        plt.plot(inp.x, sp2.scalar, 'r')
        plt.savefig(inp.filename + '2_plot.png')
        plt.clf()
        #plt.plot(inp.x, sp3.scalar, 'g')
        #plt.savefig(inp.filename + '3_plot.png')



#THIS IS TO READ CSV FILES
#            abc = (pd.read_csv(inp.filename + '_dim.csv',header=None,index_col=False)).values
#            abc.shape=(2)
#            df=(pd.read_csv(inp.filename + '_anglr.csv',header=None,index_col=False)).values
#            df.shape=(abc[0],abc[1])
