## Filename: fVUpDown.py
## Date: 22/02/24
## Author: Itelcontar
## Description: Python3 module for scalar field plotting of Y (in particular, but others are possible) velocities, above and below an impeller plane.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. cellSizeArray - Numpy array of integers corresponding to desired cell size (in mm) for each experiment in 'exp'.
# 3. component - Velocity component desired for plots. Y is recommended, and default, if not specified. Takes form 'vY','vX','vZ' or 'V' or 'VXZ'. The latter is a resultant of only X & Z.
# 4. spuriousFilter - Integer, specifying minimum number of velocity data points in a pixel for the pixel to be included in the final plot. Removes 'spurious' velocities, with insufficient data for accurate averaging, or valid assumption of ergodicity.
# 5. cutPoint - Most likely impeller height (in units of tank Diameters). E.g. 0.333.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from FGlobalSettings import globalData

### FUNCTION ###

def fVUpDown(exp, cellSizeArray, component, spuriousFilter, cutPoint):

    # Running through 'exp'
    for i in range(0,len(exp)):
        # By default, all position data is normalised to tank diameter.
        TankD = globalData[exp[i]]['TankD']
        time = globalData[exp[i]]['time']
        posX = np.divide(globalData[exp[i]]['posX'],TankD)
        posY = np.divide(globalData[exp[i]]['posY'],TankD)
        posZ = np.divide(globalData[exp[i]]['posZ'],TankD)
        name = globalData[exp[i]]['name']
        cellSize = cellSizeArray[i]

        # Selecting Velocity Component
        if component.lower() == 'vx':
            V = globalData[exp[i]]['vX']
            nameV = 'X'
        elif component.lower() == 'vz':
            V = globalData[exp[i]]['vZ']
            nameV = 'Z'
        elif component.lower() == 'v':
            V = np.sqrt(globalData[exp[i]]['vZ']**2 + globalData[exp[i]]['vY']**2 + globalData[exp[i]]['vX']**2)
            nameV = 'Magnitude'
        elif component.lower() == 'vxz':
            V = np.sqrt(globalData[exp[i]]['vZ']**2 + globalData[exp[i]]['vZ']**2)
            nameV = 'XZ Magnitude'
        else:
            V = globalData[exp[i]]['vY']
            nameV = 'Y'
            
        # Axis Labels
        xlab = 'X Position $(r/R)$'
        ylab = 'Z Position $(r/R)$'
        
        # No cells in Grid.
        nx = int(np.ceil(np.nanmax(posX)*TankD / (cellSize)))
        nz = int(np.ceil(np.nanmax(posZ)*TankD / (cellSize)))

        # Creating Grids - All same
        scalGridAbove = np.zeros((nz, nx))
        ngridAbove = np.zeros((nz, nx))
        scalGridBelow = np.zeros((nz, nx))
        ngridBelow = np.zeros((nz, nx))
        
        for i in range(len(V)):
            x, y, z, v = posX[i], posY[i], posZ[i], V[i]

            # Exclude data points with NaN values
            if not np.isnan(x) and not np.isnan(z) and not np.isnan(v):
                if y > cutPoint: 
                    ix = int(x *TankD/ (cellSize))
                    iz = int(z *TankD/ (cellSize))

                    ngridAbove[iz, ix] += 1  # Add one to the count
                    scalGridAbove[iz, ix] += v  # Adding the velocity to the velocity count                
                else:
                    ix = int(x *TankD/ (cellSize))
                    iz = int(z *TankD/ (cellSize))

                    ngridBelow[iz, ix] += 1  # Add one to the count
                    scalGridBelow[iz, ix] += v  # Adding the velocity to the velocity count                

        # To remove spurious velocities and average data.      
        scalGridAbove[ngridAbove <= spuriousFilter] = np.nan
        scalGridBelow[ngridBelow <= spuriousFilter] = np.nan
        scalGridAbove = np.nan_to_num(scalGridAbove / ngridAbove)
        scalGridBelow = np.nan_to_num(scalGridBelow / ngridBelow)

        # Use vmin and vmax to set color scale limits
        vmin = min(np.nanmin(scalGridAbove), np.nanmin(scalGridBelow))*0.9
        vmax = max(np.nanmax(scalGridAbove), np.nanmax(scalGridBelow))*0.9
        vmax = max(abs(vmax),abs(vmin))
        AR = 1

        cmap = mpl.colormaps.get_cmap('RdBu')
        cmap.set_bad(color='red')
        
        # Plotting
        fig, ax = plt.subplots(1, 2, figsize=(12,5))
        im1 = ax[0].imshow(scalGridAbove, extent=[-np.nanmax(posX), np.nanmax(posX), -np.nanmax(posZ), np.nanmax(posZ)], cmap=cmap, vmin=-vmax, vmax=vmax)
        im2 = ax[1].imshow(scalGridBelow, extent=[-np.nanmax(posX), np.nanmax(posX), -np.nanmax(posZ), np.nanmax(posZ)], cmap=cmap, vmin=-vmax, vmax=vmax)

        # Create a common colorbar for both plots
        cbar = fig.colorbar(im2, ax=ax, orientation='vertical', fraction=0.04, pad=0.025)
        cbar.set_label("mm/s", fontsize=12)
        
        # Plot Aesthetics
        ax[0].set(ylabel=ylab, xlabel=xlab)
        ax[0].set_aspect(AR)
        ax[1].set(ylabel=ylab, xlabel=xlab)
        ax[1].set_aspect(AR)
        ax[0].grid(False)
        ax[1].grid(False)
        ax[0].set_title("Above Impeller Plane", fontsize=10)
        ax[1].set_title("Below Impeller Plane", fontsize=10)
        fig.suptitle("Depth Averaged " + nameV + " Velocity Components: " + name)
