## Filename: FJValueInvestigate.py
## Date: 03/03/24
## Author: Itelcontar
## Description: Python3 module for regenerating velocities from positions, a Priori, and repeating the process with different JValues, to identify the impact on results.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. jValueRange - Integer list or np array, containing all JValues to be investigated.
# 3. plots - Boolean numpy array for generating various plots. [Velocity Distributions, Velocity Fields].

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData

### FUNCTION ###

def fJValueInvestigate(exp,jValueRange,plots):

    ### INTERNAL VARIABLES ###
    # These should not need changing, but are a good starting point for troubleshooting.
    vMax = 700 # Maximum velocity that is not ignored & removed, in mm/s.
    minRegressionPoints = 2 # Minimum number of data points required in jValue range for linear regression to take place. Recommended: 2.
    normVectors = False # Normalise vectors, or not - if normalised all vectors are equal in size. No impact on numerical values, only plot appearance.
    addTitle = True # Add title or not, to final plot.
    vWidth = 0.003 # Width of vectors in plot.
    bins = 100 # Bins for VDistrib plot.
    cellSize = 9 # Cellsize for Velocity field plot.
    spuriousFilter = 5 # FIlter to remove cells with fewer than ** velocity data points in velocity field plot.
    
    # Gathering Data
    exp = exp[0] # Only one experiment considered for this file
    name = globalData[exp]['name']
    print("Starting run " + name + " ... this may take a while...")
    np.seterr(divide='ignore',invalid='ignore') # Divide 0 error protection
    
    # Selecting Velocities from Dictionary, and normalising. 
    posX = np.divide(globalData[exp]['posX'],globalData[exp]['TankD'])
    posY = np.divide(globalData[exp]['posY'],globalData[exp]['TankD'])
    posZ = np.divide(globalData[exp]['posZ'],globalData[exp]['TankD'])
    TankD = globalData[exp]['TankD']
    time = globalData[exp]['time']
    vMax = vMax/TankD
    cellSize = cellSize
    xlab = ['Magnitude (TankD/s)','$V_{x}$ (TankD/s)','$V_{y}$ (TankD/s)','$V_{z}$ (TankD/s)'] 

    # Defining an internal function to produce velocity plots for given J Value
    def fJValuePlots(time,posX,posY,posZ,jValue,plots,count,cellSize,vMax,minRegressionPoints,vWidth):
        # Setting up empty velocity arrays
        vX = np.full(len(posX), np.NaN)
        vY = np.full(len(posY), np.NaN)
        vZ = np.full(len(posZ), np.NaN)

        # Linear regression - velocity calculations
        for n in range(jValue, len(time) - jValue):
            # Mask to avoid NaNs
            mask = ~np.isnan(posX[n - jValue:n + jValue]) & ~np.isnan(posY[n - jValue:n + jValue]) & ~np.isnan(posZ[n - jValue:n + jValue])

            # Requiring a mininum number of data points in our jValue ranges to carry out regression; 2 is recommended
            if np.sum(mask) >= minRegressionPoints:  # At least two non-NaN values required for regression
                polyX = np.polyfit(time[n - jValue:n + jValue][mask], posX[n - jValue:n + jValue][mask], 1)
                polyY = np.polyfit(time[n - jValue:n + jValue][mask], posY[n - jValue:n + jValue][mask], 1)
                polyZ = np.polyfit(time[n - jValue:n + jValue][mask], posZ[n - jValue:n + jValue][mask], 1)

                # Finding the velocity (i.e. gradient)
                vX[n] = polyX[0]
                vY[n] = polyY[0]
                vZ[n] = polyZ[0]

        # Trimming off velocities outside specified range
        vX = np.where(abs(vX) > vMax, np.nan, vX)
        vY = np.where(abs(vY) > vMax, np.nan, vY)    
        vZ = np.where(abs(vZ) > vMax, np.nan, vZ)
            
        ### PLOTTING ###
        vMag = np.sqrt(vX**2 + vY**2 + vZ**2)

        if plots[0] == True:
            # Initialising distribution figure
            fig0, ax0 = plt.subplots(4,1,figsize = (6,9), tight_layout = True)
            chartAlpha = 1
            # Plot 1 - Distributions
            ax0[0].hist(vMag, bins=bins, density=True, alpha=chartAlpha, label=name)
            ax0[1].hist(vX, bins=bins, density=True, alpha=chartAlpha, label=name)
            ax0[2].hist(vY, bins=bins, density=True, alpha=chartAlpha, label=name)
            ax0[3].hist(vZ, bins=bins, density=True, alpha=chartAlpha, label=name)

            # Plot Aesthetics
            for ax in ax0:
                ax.set_facecolor('white')
                ax.set_xlim(-1.2,1.2)
                ax.set_ylim(0,1.4)
            ax0[0].set_xlim(0,2.2)

            if addTitle == True:
                fig0.suptitle('J Value: ' + str(jValue), fontsize=14)

            ax0[0].set(xlabel=xlab[0], ylabel='Probability Density')
            ax0[1].set(xlabel=xlab[1], ylabel='Probability Density')
            ax0[2].set(xlabel=xlab[2], ylabel='Probability Density')
            ax0[3].set(xlabel=xlab[3], ylabel='Probability Density')
            ax0[0].legend()
            plt.savefig('1Jplot'+str(count)+'.png')
            
        if plots[0] == True:
            # Plot 2 - Velocity Fields (X-Y plane)
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
            ax.quiver(scaleX, scaleY, vxGrid, vyGrid, color='black', width=vWidth)

            # Plot Aesthetics
            ax.set(ylabel='$(h/H)$', xlabel='$(r/R)$')
            ax.set_aspect(AR)
            
            ax.set_facecolor('white')
            if addTitle == True:
                fig.suptitle('J Value: ' + str(jValue), fontsize=14) 

            # Saving Figure
            plt.savefig('2Jplot'+str(count)+'.png')
            
    ### ### ### ### END OF INTERNAL FUNCTION ### ### #### ### ### #### ##
    count = 0
    for i in jValueRange:
        count = count + 1
        fJValuePlots(time,posX,posY,posZ,i,plots,count,cellSize,vMax,minRegressionPoints,vWidth)
