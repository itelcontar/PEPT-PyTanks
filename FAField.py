## Filename: FAField.py
## Date: 12/03/24
## Author: Itelcontar
## Description: Python3 module for plotting acceleration fields.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. cellSize - Integer width of velocity field cells in mm.
# 3. angles - Integer numpy array of same dimensions as exp. Values either 1,2 or 3. 1 is X-Y plane, 2 is Z-Y plane, and 3 is X-Z plane. Note PEPT covention holds Y as vertical, not Z.
# 4. spuriousFilter - Integer, specifying minimum number of velocity data points in a pixel for the pixel to be included in the final plot. Removes 'spurious' accelerations, with insufficient data for accurate averaging, or valid assumption of ergodicity.
# 5. plots - Boolean array indicated which plots to create [Acceleration Distirbutions, Acceleration Fields]

### INTERNAL VARIABLES ###
# These should not need changing, but are a good starting point for troubleshooting.
aMax = 1000 # Maximum acceleration that is not ignored & removed, in mm/s/s.
minRegressionPoints = 2 # Minimum number of data points required in jValue range for linear regression to take place. Recommended: 2.
normVectors = False # Normalise vectors, or not - if normalised all vectors are equal in size. No impact on numerical values, only plot appearance.
addTitle = False # Add title or not, to final plot.
vWidth = 0.003 # Width of vecotrs in plot.
bins = 100 # Bins for ADistrib plot.
cellSize = 30 # Cellsize for Acceleration field plot.
spuriousFilter = 3 # Filter to remove cells with fewer than ** velocity data points in velocity field plot.
jValue = 10 # Integer, describing No. of datapoints either side of a position to include in linear regression for acceleration calculation.

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData
from mpl_toolkits.mplot3d import Axes3D

### FUNCTION ###

def fAField(exp, cellSize, angles, spuriousFilter,plots):

    # Running through each experiment specified
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

        # Setting up empty arrays
        aX = np.full(len(posX),np.NaN)
        aY = np.full(len(posX),np.NaN)
        aZ = np.full(len(posX),np.NaN)

        # Linear regression - acceleration calculations
        for n in range(jValue, len(time) - jValue):
            # Mask to avoid NaNs
            mask = ~np.isnan(vX[n - jValue:n + jValue]) & ~np.isnan(vY[n - jValue:n + jValue]) & ~np.isnan(vZ[n - jValue:n + jValue])

            # Requiring a mininum number of data points in our jValue ranges to carry out regression; 2 is recommended
            if np.sum(mask) >= minRegressionPoints:  # At least two non-NaN values required for regression
                polyX = np.polyfit(time[n - jValue:n + jValue][mask], vX[n - jValue:n + jValue][mask], 1)
                polyY = np.polyfit(time[n - jValue:n + jValue][mask], vY[n - jValue:n + jValue][mask], 1)
                polyZ = np.polyfit(time[n - jValue:n + jValue][mask], vZ[n - jValue:n + jValue][mask], 1)

                # Finding the acceleration (i.e. gradient)
                aX[n] = polyX[0]
                aY[n] = polyY[0]
                aZ[n] = polyZ[0]
        
        # Trimming off acceleration outside specified range
        aX = np.where(abs(aX) > aMax, np.nan, aX)
        aY = np.where(abs(aY) > aMax, np.nan, aY)    
        aZ = np.where(abs(aZ) > aMax, np.nan, aZ)

        ### PLOTTING ###
        aMag = np.sqrt(aX**2 + aY**2 + aZ**2)

        if plots[0] == True:
            # Initialising distribution figure
            fig0, ax0 = plt.subplots(4,1,figsize = (6,9), tight_layout = True)
            chartAlpha = 1
            
            # Plot 1 - Distributions
            ax0[0].hist(aMag, bins=bins, density=True, alpha=chartAlpha, label=name)
            ax0[1].hist(aX, bins=bins, density=True, alpha=chartAlpha, label=name)
            ax0[2].hist(aY, bins=bins, density=True, alpha=chartAlpha, label=name)
            ax0[3].hist(aZ, bins=bins, density=True, alpha=chartAlpha, label=name)

            # Plot Aesthetics
            for ax in ax0:
                ax.set_facecolor('white')

            if addTitle == True:
                fig0.suptitle('Acceleration Distributions : ' + name, fontsize=14)

            ax0[0].set(xlabel='Acceleration $(mm/s^{2})$', ylabel='Probability Density')
            ax0[1].set(xlabel='Acceleration $(mm/s^{2})$', ylabel='Probability Density')
            ax0[2].set(xlabel='Acceleration $(mm/s^{2})$', ylabel='Probability Density')
            ax0[3].set(xlabel='Acceleration $(mm/s^{2})$', ylabel='Probability Density')
            ax0[0].legend()
            
        if plots[1] == True:
            # Plot 2 - Acceleration Fields
                # Firstly, identifying and adjusting angles
            if angles[i] == 2:
                posX = posZ  # Z-Y, alternative side view
                aX = aZ
                xlab = '$(r/R)$'
                ylab = '$(h/H)$'
                yisy = True
            elif angles[i] == 3:
                posY = posZ  # Top down view
                aY = aZ
                xlab = '$(r/R)$'
                ylab = '$(r/R)$'
                yisy = False
            else:
                # Do nothing
                xlab = '$(r/R)$'
                ylab = '$(h/H)$'
                yisy = True

            # Numebr of cells    
            nx = int(np.ceil(np.nanmax(posX)*TankD / (cellSize)))
            ny = int(np.ceil(np.nanmax(posY)*TankD / (cellSize)))

            # Creating grids
            axGrid = np.zeros((ny, nx))
            ayGrid = np.zeros((ny, nx))
            aMagGrid = np.zeros((ny, nx))
            nGrid = np.zeros((ny, nx))

            for i in range(len(aX)):
                x, y, ax, ay, A = posX[i], posY[i], aX[i], aY[i], aMag[i]

                # Exclude data points with NaN values
                if not np.isnan(x) and not np.isnan(y) and not np.isnan(ax) and not np.isnan(ay):
                    ix = int(x*TankD/ (cellSize))
                    iy = int(y*TankD/ (cellSize))

                    nGrid[iy, ix] += 1  # Add one to the count
                    axGrid[iy, ix] += ax  # Adding the velocity to the velocity count
                    ayGrid[iy, ix] += ay
                    aMagGrid[iy, ix] += A
                
            # To remove spurious velocities
            axGrid[nGrid <= spuriousFilter] = np.nan
            ayGrid[nGrid <= spuriousFilter] = np.nan

            axGrid = np.nan_to_num(axGrid / nGrid)  # Averaging
            ayGrid = np.nan_to_num(ayGrid / nGrid)
            aMagGrid = np.nan_to_num(aMagGrid / nGrid)
            
            # Centering on 0, and normalising to tank radius. NB action depends on angle.
            ymin = 0
            ymax = 1
            yDiff = 1/ny
            AR = 2 # Aspect Ratio
            
            # Defining position grids for plotting
            scaleX, scaleY = np.meshgrid(
                np.arange(-1, 1, 2/nx),
                np.arange(ymin, ymax, yDiff))

            # Define figure
            fig, ax = plt.subplots(figsize=(8, 8), tight_layout=True)
            ax.quiver(scaleX, scaleY, axGrid, ayGrid, color='black', width=vWidth)

            # Plot Aesthetics
            ax.set(ylabel='$(h/H)$', xlabel='$(r/R)$')
            ax.set_aspect(AR)
            ax.set_facecolor('white')
            if addTitle == True:
                fig.suptitle('Acceleration Field : ' + name, fontsize=14)

            # Contour Field
            fig1, ax1 = plt.subplots(figsize=(8, 8), tight_layout=True)
            stream = ax1.streamplot(scaleX, scaleY, axGrid, ayGrid, color='black')
        
            # Plot Aesthetics
            ax1.set(ylabel='$(h/H)$', xlabel='$(r/R)$')
            ax1.set_aspect(AR)
            ax1.set_facecolor('white')

            # Scalar Field
            fig3, ax3 = plt.subplots(figsize=(8, 8), tight_layout=True)
            im = ax3.imshow(aMagGrid, origin = 'lower', extent=[np.nanmin(posX),np.nanmax(posX),ymin,ymax], cmap = 'plasma')
            ax3.set(ylabel='$(h/H)$', xlabel='$(r/R)$')
            ax3.set_aspect(AR)
            colorbar = plt.colorbar(im, ax = ax3)
            colorbar.set_label("mm/s", fontsize=12)

            ###  3D Plot ###

            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111, projection='3d')
            surf = ax2.plot_surface(scaleX, scaleY, np.sqrt(axGrid**2 + ayGrid**2), cmap='jet')
