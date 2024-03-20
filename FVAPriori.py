## Filename: FVAPriori.py
## Date: 03/03/24
## Author: Itelcontar
## Description: Python3 module for regenerating velocities from positions, a Priori, and overwriting the source file with the new velocities.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. jValue - Integer, describing No. of datapoints either side of a position to include in linear regression for velocity calculation.
# 3. plots - Boolean numpy array for generating various plots. [Velocity Distributions, Velocity Fields].

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData

### FUNCTION ###

def fVAPriori(exp,jValue,plots):

    ### INTERNAL VARIABLES ###
    # These should not need changing, but are a good starting point for troubleshooting.
    vMax = 700 # Maximum velocity that is not ignored & removed, in mm/s.
    minRegressionPoints = 2 # Minimum number of data points required in jValue range for linear regression to take place. Recommended: 2.
    normVectors = False # Normalise vectors, or not - if normalised all vectors are equal in size. No impact on numerical values, only plot appearance.
    addTitle = True # Add title or not, to final plot.
    vWidth = 0.003 # Width of vecotrs in plot.
    bins = 100 # Bins for VDistrib plot.
    cellSize = 9 # Cellsize for Velocity field plot.
    spuriousFilter = 3 # FIlter to remove cells with fewer than ** velocity data points in velocity field plot.
    writeData = False # Whether to write the calculated velocity data back to a file. 
    fileName = "" # The filename to write the data to. Set to the same file as the exp data.

    # Initialising distribution figure
    fig0, ax0 = plt.subplots(4,1,figsize = (7,9), tight_layout = True)
    chartAlpha = 1
    
    # Running through each experiment specified
    for i in exp:

        # Collecting data

        # Selecting Velocities from Dictionary, and normalising. 
        posX = np.divide(globalData[i]['posX'],globalData[i]['TankD'])
        posY = np.divide(globalData[i]['posY'],globalData[i]['TankD'])
        posZ = np.divide(globalData[i]['posZ'],globalData[i]['TankD'])
        TankD = globalData[i]['TankD']
        time = globalData[i]['time']
        vMax = vMax/TankD
        xlab = ['Magnitude (TankD/s)','$V_{x}$ (TankD/s)','$V_{y}$ (TankD/s)','$V_{z}$ (TankD/s)'] 
        cellSize = cellSize/TankD
            
        name = globalData[i]['name']
        print("Starting run " + name + " ... this may take a while...")
        np.seterr(divide='ignore',invalid='ignore') # Divide 0 error protection

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
            # Plot 1 - Distributions
            ax0[0].hist(vMag, bins=bins, density=True, alpha=chartAlpha, label=name)
            ax0[1].hist(vX, bins=bins, density=True, alpha=chartAlpha, label=name)
            ax0[2].hist(vY, bins=bins, density=True, alpha=chartAlpha, label=name)
            ax0[3].hist(vZ, bins=bins, density=True, alpha=chartAlpha, label=name)

            # Plot Aesthetics
            for ax in ax0:
                ax.set_facecolor('white')
                ax.set_xlim(-1.2,1.2)
                ax.set_ylim(0,1.5)
            ax0[0].set_xlim(0,2.2)

            if addTitle == True:
                fig0.suptitle("Velocity Distributions"  + name + ', J Value: ' + str(jValue), fontsize = 14)

            ax0[0].set(xlabel=xlab[0], ylabel='Probability Density')
            ax0[1].set(xlabel=xlab[1], ylabel='Probability Density')
            ax0[2].set(xlabel=xlab[2], ylabel='Probability Density')
            ax0[3].set(xlabel=xlab[3], ylabel='Probability Density')
            ax0[0].legend()
            
            # Saving Figure
            plt.savefig(name+"1_J_"+str(jValue)+'.png')
            
        if plots[0] == True:
            # Plot 2 - Velocity Fields (X-Y plane)
            # No cells in grid
            nx = int(np.ceil(np.nanmax(posX) / (cellSize)))
            ny = int(np.ceil(np.nanmax(posY) / (cellSize)))
            print(nx)
            print(ny)

            # Creating grids
            vxGrid = np.zeros((ny, nx))
            vyGrid = np.zeros((ny, nx))
            nGrid = np.zeros((ny, nx))

            for i in range(len(vX)):
                x, y, vx, vy = posX[i], posY[i], vX[i], vY[i]

                # Exclude data points with NaN values
                if not np.isnan(x) and not np.isnan(y) and not np.isnan(vx) and not np.isnan(vy):
                    ix = int(x / (cellSize))
                    iy = int(y / (cellSize))

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

            # Normaise the plot vectors?
            if normVectors == True:
                normGrid = np.sqrt(vxGrid ** 2 + vyGrid ** 2)
            else:
                normGrid = 1

            # Define figure
            fig, ax = plt.subplots(figsize=(8, 8), tight_layout=True)
            ax.quiver(scaleX, scaleY, vxGrid / normGrid, vyGrid / normGrid, color='black', width=vWidth)

            # Plot Aesthetics
            ax.set(ylabel='$(h/H)$', xlabel='$(r/R)$')
            ax.set_aspect(AR)
            
            ax.set_facecolor('white')
            if addTitle == True:
                fig.suptitle("Velocity Field: "  + name + ', J Value: ' + str(jValue), fontsize=14) 

            plt.savefig(name+"2_J_"+str(jValue)+'.png')
