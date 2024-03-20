## Filename: FVDistrib.py
## Date: 22/12/23
## Author: Itelcontar
## Description: Python3 module for creating velocity distribution plots.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. bins - Number of bins to include in histogram plots.
# 3. norm - Determines whether to 'normalise' velocities. Options are 'TankD' (Tank Diameter), 'TipSpeed' (Impeller Tip Velocity), or nothing (any other value).
        # Not case-sensitive. Both diameter, and tip speed calculated from data stored in globalData[].
# 4. plots - Boolean numpy array for generating various plots. [Histogram, Line Chart].     

### INTERNAL VARIABLES ###
# These should not need changing, but are a good starting point for troubleshooting.
XaxisCompress = 0.9 # Cuts data to 0.** of the maximum visible data. Improves appearance of gausians. Change if plot view is missing data.

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData

### FUNCTION ###

def fVDistrib(exp, bins, norm, plots):
    
    # Initialising figures
    fig1, ax1 = plt.subplots(4,1,figsize = (6,9), tight_layout = True)
    fig2, ax2 = plt.subplots(4,1,figsize = (6,9), tight_layout = True)
    chartAlpha = 1/len(exp)
    
    # Running through 'exp'
    for i in exp:
        if norm.lower() == 'tankd':
            # Selecting velocities from Dictionary, and normalising 
            vX = np.divide(globalData[i]['vX'],globalData[i]['TankD']) 
            vY = np.divide(globalData[i]['vY'],globalData[i]['TankD'])
            vZ = np.divide(globalData[i]['vZ'],globalData[i]['TankD'])
            RPM = globalData[i]['RPM']
            title = "Tank Velocity Probability Density Distributions - Normalised with Tank Diameter"
            xlab = ['Magnitude (TankD/s)','$V_{x}$ (TankD/s)','$V_{y}$ (TankD/s)','$V_{z}$ (TankD/s)'] 

        elif norm.lower() == 'tipspeed':
            # Finding impeller tip speed
            radius = globalData[i]['TankD'] / 2
            RPM = globalData[i]['RPM']
            tipSpeed = (1/3) * radius * 2 * np.pi * RPM / 60 # (1/3) * 2piR * RPM / 60 => mm/s tip speed            

            # Selecting velocities from Dictionary, and normalising
            vX = np.divide(globalData[i]['vX'],tipSpeed) 
            vY = np.divide(globalData[i]['vY'],tipSpeed)
            vZ = np.divide(globalData[i]['vZ'],tipSpeed)
            title = "Tank Velocity Probability Density Distributions - Normalised with Impeller Tip Velocity"
            xlab = ['Magnitude ($v/V_{Tip}$)','$v_{x}/V_{Tip}$','$v_{y}/V_{Tip}$','$v_{z}/V_{Tip}$'] 

        else:
            # Selecting velocities from dictionary 
            vX = globalData[i]['vX']
            vY = globalData[i]['vY']
            vZ = globalData[i]['vZ']
            RPM = globalData[i]['RPM']
            title = "Tank Velocity Probability Density Distributions"
            xlab = ['Magnitude (mm/s)','$V_{x}$ (mm/s)','$V_{y}$ (mm/s)','$V_{z}$ (mm/s)'] 
            
        vMag = np.sqrt(vX**2 + vY**2 + vZ**2)
        iname = globalData[i]['name']


        ## PLOT 1 ##
        if plots[0] == True:
            # Probability distribution
            ax1[0].hist(vMag, bins=bins, density=True, alpha=chartAlpha, label=iname + ", " + str(RPM) + " RPM")
            ax1[1].hist(vX, bins=bins, density=True, alpha=chartAlpha, label=iname + ", " + str(RPM) + " RPM")
            ax1[2].hist(vY, bins=bins, density=True, alpha=chartAlpha, label=iname + ", " + str(RPM) + " RPM")
            ax1[3].hist(vZ, bins=bins, density=True, alpha=chartAlpha, label=iname + ", " + str(RPM) + " RPM")

            # Compressing X axis to improve view. Disable if encountering viewing issues
            currentXLimMax = max(ax1[1].get_xlim()[1], ax1[2].get_xlim()[1], ax1[3].get_xlim()[1])
            currentXLimMin = min(ax1[1].get_xlim()[0], ax1[2].get_xlim()[0], ax1[3].get_xlim()[0])
            currentXLim = max(abs(currentXLimMin),abs(currentXLimMax))
            newXLimMax = XaxisCompress * currentXLim
            ax1[1].set_xlim(-newXLimMax,newXLimMax)
            ax1[2].set_xlim(-newXLimMax,newXLimMax)
            ax1[3].set_xlim(-newXLimMax,newXLimMax)
            ax1[0].set_xlim(0,XaxisCompress*ax1[0].get_xlim()[1])

            # Plot Aesthetics
            for ax in ax1:
                ax.set_facecolor('white')
            ax1[0].set(xlabel=xlab[0], ylabel='P()')
            ax1[1].set(xlabel=xlab[1], ylabel='P()')
            ax1[2].set(xlabel=xlab[2], ylabel='P()')
            ax1[3].set(xlabel=xlab[3], ylabel='P()')
            ax1[0].legend()
            fig1.suptitle(title, fontsize=10)

        ## PLOT 2 ##
        if plots[1] == True:
            # Lines chart, not Bars!
            nMag, binsMag, _ = ax2[0].hist(vMag, bins=bins, density=True, alpha=0)
            nx, binsX, _ = ax2[1].hist(vX, bins=bins, density=True, alpha=0)
            ny, binsY, _ = ax2[2].hist(vY, bins=bins, density=True, alpha=0)
            nz, binsZ, _ = ax2[3].hist(vZ, bins=bins, density=True, alpha=0)
        
            # Calculate bin centers
            binXCenters = (binsX[:-1] + binsX[1:]) / 2
            binYCenters = (binsY[:-1] + binsY[1:]) / 2
            binZCenters = (binsZ[:-1] + binsZ[1:]) / 2
            binMagCenters = (binsMag[:-1] + binsMag[1:]) / 2

            # Plot the line using bin centers and histogram values
            ax2[0].plot(binMagCenters, nMag, label=iname + ", " + str(RPM) + " RPM")
            ax2[1].plot(binXCenters, nx, label=iname + ", " + str(RPM) + " RPM")
            ax2[2].plot(binYCenters, ny, label=iname + ", " + str(RPM) + " RPM")
            ax2[3].plot(binZCenters, nz, label=iname + ", " + str(RPM) + " RPM")

            # Compressing X axis to improve view. Disable if encountering viewing issues
            currentXLimMax = max(ax2[1].get_xlim()[1], ax2[2].get_xlim()[1], ax2[3].get_xlim()[1])
            currentXLimMin = min(ax2[1].get_xlim()[0], ax2[2].get_xlim()[0], ax2[3].get_xlim()[0])
            currentXLim = max(abs(currentXLimMin),abs(currentXLimMax))
            newXLimMax = XaxisCompress * currentXLim
            ax2[1].set_xlim(-newXLimMax,newXLimMax)
            ax2[2].set_xlim(-newXLimMax,newXLimMax)
            ax2[3].set_xlim(-newXLimMax,newXLimMax)
            ax2[0].set_xlim(0,XaxisCompress*ax2[0].get_xlim()[1])
            
            # Plot Aesthetics
            for ax in ax2:
                ax.set_facecolor('white')
            ax2[0].set(xlabel=xlab[0], ylabel='P()')
            ax2[1].set(xlabel=xlab[1], ylabel='P()')
            ax2[2].set(xlabel=xlab[2], ylabel='P()')
            ax2[3].set(xlabel=xlab[3], ylabel='P()')
            ax2[0].legend()
            fig2.suptitle(title, fontsize=10)
