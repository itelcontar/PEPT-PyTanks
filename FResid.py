## Filename: fResid.py
## Date: 22/12/23
## Author: Itelcontar
## Description: Python3 module for creating Residency Distributions.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. cellSize - Integer width of residency field cells in mm.
# 3. angles - Integer numpy array of same dimensions as exp. Values either 1,2 or 3. 1 is X-Y plane, 2 is Z-Y plane, and 3 is X-Z plane. Note PEPT covention holds Y as vertical, not Z.

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData

### FUNCTION ### 

def fResid(exp,cellSizeArray,angleArray):

    for n in range(0,len(exp)):
        # Data gathering
        posX = globalData[exp[n]]['posX']
        posY = globalData[exp[n]]['posY']
        posZ = globalData[exp[n]]['posZ']
        time = globalData[exp[n]]['time']
        name = globalData[exp[n]]['name']
        cellSize = cellSizeArray[n]
        angle = angleArray[n]
        
        # First check and change angle:
        if angle == 2:
            posX = posZ # Alternative side view
            xlab = 'Z Position $(mm)$'
            ylab = 'Y Position $(mm)$'
            ignoreCentering = False
        elif angle == 3:
            posY = posZ # Top down view
            xlab = 'X Position $(mm)$'
            ylab = 'Z Position $(mm)$'
            ignoreCentering = True
        else:
            # Do nothing
            xlab = 'X Position $(mm)$'
            ylab = 'Y Position $(mm)$'
            ignoreCentering = False
            
        posMid = np.zeros((len(posX) - 1, 3))
        dt = np.zeros(len(posX) - 1)

        for i in range(0,len(posX)-1): # finding average positions to clean data a bit
            posMid[i,0] = (posX[i] + posX[i+1]) / 2 
            posMid[i,1] = (posY[i] + posY[i+1]) / 2
            posMid[i,2] = (posZ[i] + posZ[i+1]) / 2
            dt[i] = time[i+1] - time[i]

        # No cells in Grid.
        nx = int(np.ceil(np.nanmax(posX) / cellSize))
        ny = int(np.ceil(np.nanmax(posY) / cellSize))

        # Creating Grid
        Resgrid = np.zeros((ny,nx))

        # Mask Check, to avoid NaNs
        ResMask = np.logical_not(np.logical_or( np.isnan(posX), np.isnan(posY)))

        # Loop through and add up time spent in each cell
        for i in range(0,len(posX) - 1):
            if ResMask[i] == True: # NB true means a clean data point with no NaNs.
                ix = int(posX[i] / cellSize)
                iy = int(posY[i] / cellSize)
                Resgrid[iy,ix] = Resgrid[iy,ix] + dt[i]  

        norm_on = True
        if norm_on == True:
            norm = Resgrid / Resgrid.max()
        else:
            norm = 1

        # Plotting
        fig4, ax4 = plt.subplots(figsize = (8,8), tight_layout = True)
        im = ax4.imshow(Resgrid,origin='lower',extent=[0,np.nanmax(posX),0,np.nanmax(posY)], cmap = 'jet')

        # Plot Aesthetics
        ax4.set(ylabel = ylab, xlabel = xlab)
        ax4.set_aspect(1)
        ax4.grid(False)
        fig4.suptitle("Residency Distribution: " + name, fontsize = 14)
        colorbar = plt.colorbar(im, ax = ax4)

