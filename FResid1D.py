## Filename: FResid1D.py
## Date: 15/01/24
## Author: Itelcontar
## Description: Python3 module for creating Residency Distribution in 1 dimension.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. cellSize - Integer width of residency field cells in mm.
# 3. angles - Which axis should be used for the 1D plot? String: 'posX', 'posY' or 'posZ'.
# 4. norm - Boolean True/False. Indicates if plots should be 'normalised' to tank Diameter. 
# 5. nPos - Number of points to use: more, greater noise, accuracy and time; less lower granularity & faster.

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData

### FUNCTION ###

def fResid1D(exp, angle, norm, nPos):
    # Initializing figure variables
    fig1, ax1 = plt.subplots(1, 1,figsize = (7,9), tight_layout = True)    
    axMax = 0
    
    # Running through each experiment specified
    for i in exp:
        # Sorting out angle
        dim = angle[3]
        if dim == "Y":
            ignoreCentering = True
        else:
            ignoreCentering = False
        
        if norm == True:
            # Selecting position vector from Dictionary, and normalizing. 
            posN = np.divide(globalData[i][angle], globalData[i]['TankD'])
            title = "Tank Residency Distributions - Normalized with Tank Diameter"
            xlab = 'Residency Probability'
            ylab = dim + ' Position (TankD)'
        else:
            # Selecting position vector from Dictionary 
            posN = globalData[i][angle]
            title = "Tank Residency Distributions - Raw"
            xlab = 'Residency Probability'
            ylab = dim + ' Position (mm)' 

        # Collecting experiment name    
        iname = globalData[i]['name']
        
        # Now average the position
        time = globalData[i]['time']
        posMid = np.zeros((len(posN) - 1))
        dt = np.zeros(len(posN) - 1)
        for n in range(0, len(posN) - 1): # finding average positions to clean data a bit
            posMid[n] = (posN[n] + posN[n+1]) / 2 
            dt[n] = time[n+1] - time[n]

        # Creating Grid
        Resgrid = np.zeros(nPos)
        nPosMat = np.linspace(np.nanmin(posMid), np.nanmax(posMid), num = nPos)
        cellSize = nPosMat[1] - nPosMat[0] 
        
        # Mask Check, to avoid NaNs
        ResMask = np.logical_not(np.isnan(posMid))
        
        # Loop through and add up time spent in each cell
        for n in range(0, len(posMid) - 1):
            if ResMask[n] == True: # NB true means a clean data point with no NaNs.
                nx = int((posMid[n] - nPosMat[0] ) / cellSize)
                Resgrid[nx] = Resgrid[nx] + dt[n]  

        norm_on = True # HANDLE - Do we make this a probability density distribution?
        if norm_on == True:
            # Normalize the area between the plot and the y-axis
            area = np.trapz(Resgrid, np.abs(nPosMat))
            Resgrid = Resgrid / area
        
        # Correcting 0 error
        Resgrid = np.insert(Resgrid,0,0)
        nPosMat = np.insert(nPosMat,0,0)

        # Centering
        if ignoreCentering == False:
            nPosMat = nPosMat - np.nanmean(posN)
            
        # Now We Plot
        ax1.plot(Resgrid, nPosMat, label=iname)

        if np.nanmax(Resgrid) > axMax:
            axMax = np.nanmax(Resgrid) 
        
    # Plot Aesthetics
    ax1.set(ylabel=ylab, xlabel=xlab)
    fig1.suptitle(title, fontsize=14)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.set_facecolor('#f0f0f0')  # Light grey background
    ax1.legend()

    # Setting axis limits
    ax1.set_xlim([0, axMax])
