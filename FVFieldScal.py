## Filename: VFieldScal.py
## Date: 22/02/24
## Author: Itelcontar
## Description: Python3 module for creating Velocity Scalar Field Plotting.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. cellSize - Integer width of velocity field cells in mm.
# 3. angles - Integer numpy array of same dimensions as exp. Values either 1,2 or 3. 1 is X-Y plane, 2 is Z-Y plane, and 3 is X-Z plane. Note PEPT covention holds Y as vertical, not Z.
# 4. spuriousFilter - Integer, specifying minimum number of velocity data points in a pixel for the pixel to be included in the final plot. Removes 'spurious' velocities, with insufficient data for accurate averaging, or valid assumption of ergodicity.

### INTERNAL VARIABLES ###
# These should not need changing, but are a good starting point for troubleshooting.
norm_on = False  # Divide all cells by maximum cell value?

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData

### FUNCTION ###

def fVFieldScal(exp, cellSize, angles, spuriousFilter):
    # For loop for i in exp.
    for i in range(0,len(exp)):

        # Data gathering
        TankD = globalData[exp[i]]['TankD']
        time = globalData[exp[i]]['time']
        posX = np.divide(globalData[exp[i]]['posX'],TankD)
        posY = np.divide(globalData[exp[i]]['posY'],TankD)
        posZ = np.divide(globalData[exp[i]]['posZ'],TankD)
        vX = globalData[exp[i]]['vX']
        vY = globalData[exp[i]]['vY']
        vZ = globalData[exp[i]]['vZ']
        name = globalData[exp[i]]['name']

        vMag = np.sqrt(vX**2 + vY**2 + vZ**2)
        
        # Check and change angle
        if angles[i] == 2:
            posX = posZ  # Alternative side view
            vX = vZ
            xlab = '$(r/R)$'
            ylab = '$(h/H)$'
            yisy = True
        elif angles[i] == 3:
            posY = posZ  # Top down view
            vY = vZ
            xlab = '$(r/R)$'
            ylab = '$(r/R)$'
            yisy = False
        else:
            # Do nothing
            xlab = '$(r/R)$'
            ylab = '$(h/H)$'
            yisy = True

        # No cells in Grid.
        nx = int(np.ceil(np.nanmax(posX)*TankD / (cellSize)))
        ny = int(np.ceil(np.nanmax(posY)*TankD / (cellSize)))

        # Creating Grids - All same
        scalGrid = np.zeros((ny, nx))
        ngrid = np.zeros((ny, nx))

        for i in range(len(vX)):
            x, y, v = posX[i], posY[i], vMag[i]

            # Exclude data points with NaN values
            if not np.isnan(x) and not np.isnan(y) and not np.isnan(v):
                ix = int(x *TankD/ (cellSize))
                iy = int(y *TankD/ (cellSize))

                ngrid[iy, ix] += 1  # Add one to the count
                scalGrid[iy, ix] += v  # Adding the velocity to the velocity count                

        # Remove spurious velocities
        scalGrid[ngrid <= spuriousFilter] = np.nan
        scalGrid = np.nan_to_num(scalGrid / ngrid)  # Averaging

        # Before calculating position grids, understand axes limits, and centre on 0
        if yisy == True:
            posX = posX - np.nanmean(posX)
            ymax = 1
            ymin = 1 - (np.nanmax(posY)-np.nanmin(posY)) # Assumes surface has been imaged
            ydiff = 1/ny
        else:
            posX = posX - np.nanmean(posX)
            posY = posY - np.nanmean(posY)
            ymin = np.nanmin(posY)
            ymax = np.nanmax(posY)
            ydiff = 2/ny
            
        if norm_on == True:
            scalGrid = scalGrid / scalGrid.max()
                
        # Plotting
        fig, ax = plt.subplots(figsize=(8, 8), tight_layout=True)
        im = ax.imshow(scalGrid, origin = 'lower', extent=[np.nanmin(posX),np.nanmax(posX),ymin,ymax], cmap = 'plasma')

        # Aspect
        rows, cols = scalGrid.shape
        AR = cols/rows

        # Plot Aesthetics
        ax.set(ylabel=ylab, xlabel=xlab)
        ax.set_aspect(AR)
        ax.grid(False)
        fig.suptitle("Scalar Field: " + name, fontsize=14)
        colorbar = plt.colorbar(im, ax = ax)
        colorbar.set_label("mm/s", fontsize=12)

