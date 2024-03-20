## Filename: FVDistribSingle.py
## Date: 22/12/23
## Author: Itelcontar
## Description: Python3 module for creating 1x1 velocity distribution plots.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. bins - Number of bins to include in histogram plots.
# 3. norm - Determines whether to 'normalise' velocities. Options are 'TankD' (Tank Diameter), 'TipSpeed' (Impeller Tip Velocity), or nothing (any other value).
        # Not case-sensitive. Both diameter, and tip speed calculated from data stored in globalData[].
# 4. component - Determines plot velocity component. May be 'V','vZ','vY','vX'
# 5. plots - Boolean numpy array for generating various plots. [Histogram, Line Chart].     

### INTERNAL VARIABLES ###
# These should not need changing, but are a good starting point for troubleshooting.
XaxisCompress = 0.9 # Cuts data to 0.** of the maximum visible data. Improves appearance of gausians. Change if plot view is missing data.

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData

### FUNCTION ###

def fVDistribSingle(exp, bins, norm, component, plots):
    
    # Initialising figures
    lineStyles = np.array(['None','-', '--', '-.', ':', '.', ',', 'o', 'v', '^', '<', '>', 's', 'p', '+', 'x', 'D', 'd'])
    count = 0
    chartAlpha = 1/len(exp)
    if plots[0] == True:
        fig1, ax1 = plt.subplots(figsize = (6,6), tight_layout = True)
    if plots[1] == True:
        fig2, ax2 = plt.subplots(figsize = (6,6), tight_layout = True)
    
    # Running through 'exp'
    for i in exp:
        if norm.lower() ==  'tipspeed':
            # Finding impeller tip speed
            radius = globalData[i]['TankD'] / 2
            RPM = globalData[i]['RPM']
            tipSpeed = (1/3) * radius * 2 * np.pi * RPM / 60 # (1/3) * 2piR * RPM / 60 => mm/s tip speed            
            divide = tipSpeed
            xlabV = ['Magnitude ($v/V_{Tip}$)','$v_{x}/V_{Tip}$','$v_{y}/V_{Tip}$','$v_{z}/V_{Tip}$'] 
            print("Normalising to Tip Speed")
        elif norm.lower() ==  'tankd':
            # Finding impeller tip speed
            TankD = globalData[i]['TankD']
            RPM = globalData[i]['RPM']
            divide = TankD
            xlabV = ['Magnitude ($T/s$)','$v_{x}$ (T/s)','$v_{y} (T/s)$','$v_{z}$ (T/s)'] 
            print("Normalising to Tank D")
        else:
            divide = 1
            xlabV = ['Magnitude (mm/s)','$V_{x}$ (mm/s)','$V_{y}$ (mm/s)','$V_{z}$ (mm/s)'] 
        count = count + 1
        
        # Collecting Data
        iname = globalData[i]['name']
        if component.lower() == 'v':
            vX = np.divide(globalData[i]['vX'],divide) 
            vY = np.divide(globalData[i]['vY'],divide)
            vZ = np.divide(globalData[i]['vZ'],divide)
            RPM = globalData[i]['RPM']
            title = "Velocity Magnitudes"
            V = np.sqrt(vX**2 + vY**2 + vZ**2)
            xlab = xlabV[0]
        elif component.lower() == 'vx':
            V = np.divide(globalData[i]['vX'],divide)
            RPM = globalData[i]['RPM']
            title = "X Velocity Components"
            xlab = xlabV[1]
        elif component.lower() == 'vy':
            V = np.divide(globalData[i]['vY'],divide)
            RPM = globalData[i]['RPM']
            title = "Y Velocity Components"
            xlab = xlabV[2]
        elif component.lower() == 'vz':
            V = np.divide(globalData[i]['vZ'],divide)
            RPM = globalData[i]['RPM']
            title = "Z Velocity Components"
            xlab = xlabV[3]

        ## PLOT 1 ##
        if plots[0] == True:
            # Probability distribution
            ax1.hist(V, bins=bins, density=True, alpha=chartAlpha, label=iname + ", " + str(RPM) + " RPM")

            # Plot Aesthetics
            ax1.set_facecolor('white')
            ax1.set(xlabel=xlab, ylabel='P()')
            ax1.legend()
            fig1.suptitle(title, fontsize=10)

        ## PLOT 2 ##
        if plots[1] == True:
            # Lines chart, not Bars!
            nV, binsV, _ = ax2.hist(V, bins=bins, density=True, alpha=0)
            
            # Calculate bin centers
            binVCenters = (binsV[:-1] + binsV[1:]) / 2

            # Use a colormap (e.g., viridis)
            #cmap = plt.get_cmap('turbo')
            color = 'black' #cmap(0.5*count / len(exp))

            # Plot the line using bin centers and histogram values
            ax2.plot(binVCenters, nV, color=color, linestyle = lineStyles[count],label=iname + ", " + str(RPM) + " RPM")
            
            # Plot Aesthetics
            ax2.set_facecolor('white')
            #ax2.set(xlabel=xlab, ylabel='Probability Density')
            ax2.set_xlabel(xlab, fontsize=14)
            ax2.set_ylabel('Probability Density', fontsize=14)
            ax2.tick_params(axis='x', labelsize=12)  # Set font size for x-axis tick labels
            ax2.tick_params(axis='y', labelsize=12)
            
            #ax2.set_xlim(-300,300)
            legend = ax2.legend(fontsize=10.5)
            legend.get_frame().set_edgecolor('black')
            fig2.suptitle(title, fontsize=10)
