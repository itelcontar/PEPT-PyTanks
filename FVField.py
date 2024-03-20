## Filename: FVField.py
## Date: 22/12/23
## Author: Itelcontar
## Description: Python3 module for plotting velocity fields.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. cellSize - Integer width of velocity field cells in mm.
# 3. angles - Integer numpy array of same dimensions as exp. Values either 1,2 or 3. 1 is X-Y plane, 2 is Z-Y plane, and 3 is X-Z plane. Note PEPT covention holds Y as vertical, not Z.
# 4. spuriousFilter - Integer, specifying minimum number of velocity data points in a pixel for the pixel to be included in the final plot. Removes 'spurious' velocities, with insufficient data for accurate averaging, or valid assumption of ergodicity.

### INTERNAL VARIABLES ###
# These should not need changing, but are a good starting point for troubleshooting.
normVectors = False # Normalise vectors, or not - if normalised all vectors are equal in size. No impact on numerical values, only plot appearance.
addTitle = False # Add title or not, to final plot.
vWidth = 0.003 # Width of vecotrs in plot.

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData

### FUNCTION ###

def fVField(exp, cellSize, angles, spuriousFilter):

    # Running through each experiment specified
    for i in range(0,len(exp)):

        # Data gathering
        TankD = globalData[exp[i]]['TankD']
        time = globalData[exp[i]]['time']
        posX = np.divide(globalData[exp[i]]['posX'],TankD)
        posY = np.divide(globalData[exp[i]]['posY'],TankD)
        posZ = np.divide(globalData[exp[i]]['posZ'],TankD)
        vX = np.divide(globalData[exp[i]]['vX'],TankD)
        vY = np.divide(globalData[exp[i]]['vY'],TankD)
        vZ = np.divide(globalData[exp[i]]['vZ'],TankD)
        name = globalData[exp[i]]['name']
        
        # Firstly, identifying and adjusting angles
        if angles[i] == 2:
            posX = posZ  # Z-Y, alternative side view
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

        # No cells in grid
        nx = int(np.ceil(np.nanmax(posX)*TankD / (cellSize)))
        ny = int(np.ceil(np.nanmax(posY)*TankD / (cellSize)))

        # Creating grids
        vxGrid = np.zeros((ny, nx))
        vyGrid = np.zeros((ny, nx))
        nGrid = np.zeros((ny, nx))

        for i in range(len(vX)):
            x, y, vx, vy = posX[i], posY[i], vX[i], vY[i]

            # Exclude data points with NaN values
            if not np.isnan(x) and not np.isnan(y) and not np.isnan(vx) and not np.isnan(vy):
                ix = int(x *TankD/ (cellSize))
                iy = int(y *TankD/ (cellSize))

                nGrid[iy, ix] += 1  # Add one to the count
                vxGrid[iy, ix] += vx  # Adding the velocity to the velocity count
                vyGrid[iy, ix] += vy

        # To remove spurious velocities
        vxGrid[nGrid <= spuriousFilter] = np.nan
        vyGrid[nGrid <= spuriousFilter] = np.nan

        vxGrid = np.nan_to_num(vxGrid / nGrid)  # Averaging
        vyGrid = np.nan_to_num(vyGrid / nGrid)

        # Centering on 0. NB action depends on angle.
        if yisy == True:
            posX = posX - np.nanmean(posX)
            ymin = np.nanmin(posY) #0 #  # Assumes surface has been imaged
            ymax = np.nanmax(posY)
            ydiff = (ymax-ymin) / ny
            AR = 2
        else:
            posX = posX - np.nanmean(posX)
            posY = posY - np.nanmean(posY)
            ymin = np.nanmin(posY)
            ymax = np.nanmax(posY)
            ydiff = (np.nanmax(posY)  - np.nanmin(posY)  )  /ny #2/ny
            AR = 1
            
        # Defining position grids for plotting
        scaleX, scaleY = np.meshgrid(
            np.linspace(np.nanmin(posX), np.nanmax(posX), nx),
            np.linspace(ymin, ymax, ny)
        )
                
        # Normaise the plot vectors?
        if normVectors == True:
            normGrid = np.sqrt(vxGrid ** 2 + vyGrid ** 2)
        else:
            normGrid = 1
            
        # Define figure
        fig, ax = plt.subplots(figsize=(6, 6), tight_layout=True)
        ax.quiver(scaleX, scaleY, vxGrid / normGrid, vyGrid / normGrid, color='black', width=vWidth)

        # Plot Aesthetics
        rows, cols = nGrid.shape
        AR = cols/rows
        ax.set(ylabel=ylab, xlabel=xlab)
        ax.set_aspect(AR)
        ax.grid(False)

        # Adjustments
        ax.set_xlabel(xlab, fontsize=16)
        ax.set_ylabel(ylab, fontsize=16)
        ax.tick_params(axis='x', labelsize=15)  # Set font size for x-axis tick labels
        ax.tick_params(axis='y', labelsize=15)
        
        ax.set_facecolor('white')
        if addTitle == True:
            fig.suptitle("Velocity Field: " + name, fontsize=14) 
