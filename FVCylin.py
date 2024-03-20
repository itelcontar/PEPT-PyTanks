## Filename: FVCylin.py
## Date: 26/02/24
## Author: Itelcontar
## Description: Python3 module for creating cylindrical Velocities.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. jvalue - No. datapoints before and after a position to use in linear regression to find velocities.
# 3. bins - No. bins for velocity distribution plot.
# 4. norm - Normalise data to 'tipspeed' or 'tankd'.
# 5. nCellsWidth - Number of cells to divide tank radius into for velocity field plot.

### INTERNAL VARIABLES ###
# These should not need changing, but are a good starting point for troubleshooting.
normVectors = False # Normalise vectors, or not - if normalised all vectors are equal in size. No impact on numerical values, only plot appearance.
addTitle = True # Add title or not, to final plot.
vWidth = 0.005 # Width of vecotrs in plot.
spuriousFilter = 5 # Integer, specifying minimum number of velocity data points in a pixel for the pixel to be included in the final plot. Removes 'spurious' velocities, with insufficient data for accurate averaging, or valid assumption of ergodicity.

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData

### FUNCTION ###

def fVCylin(exp,jValue,bins,norm,nCellsWidth):

    # Initialising distribution figure
    fig0, ax0 = plt.subplots(4,1,figsize = (7,9), tight_layout = True)
    chartAlpha = 1/len(exp)
    # Running through each experiment specified
    for i in exp:
        
        if norm.lower() == 'tankd':
            # Selecting Velocities from Dictionary, and normalising. 
            posX = np.divide(globalData[i]['posX'],globalData[i]['TankD'])
            posY = np.divide(globalData[i]['posY'],globalData[i]['TankD'])
            posZ = np.divide(globalData[i]['posZ'],globalData[i]['TankD'])
            TankD = globalData[i]['TankD']
            time = globalData[i]['time']
            vMax = 2
            xlab = ['Vr Component (TankD/s)','Vy Component (TankD/s)','Vtheta Component (rad/s)','Velocity (TankD/s)'] 
        elif norm.lower() == 'tipspeed':
            # Selecting Velocities from Dictionary, and normalising.
            tipSpeed = (1/3) * globalData[i]['TankD'] * np.pi * globalData[i]['RPM'] / 60 # mm/s tip velocity
            RPM = globalData[i]['RPM']
            posX = globalData[i]['posX']
            posY = globalData[i]['posY']
            posZ = globalData[i]['posZ']
            TankD = globalData[i]['TankD']
            time = globalData[i]['time']
            vMax = 0.5
            xlab = ['$V_{r}/V_{Tip}$','$V_{y}/V_{Tip}$','$V_{theta}/V_{Tip theta}$','$V/V_{Tip}$']
        else:
            # Selecting Velocities from Dictionary 
            posX = globalData[i]['posX']
            posY = globalData[i]['posY']
            posZ = globalData[i]['posZ']
            TankD = globalData[i]['TankD']
            time = globalData[i]['time']
            vMax = 500* TankD/450
            #title = "Tank Velocity Probability Density Distributions"
            xlab = ['Vr Component (mm/s)','Vy Component (mm/s)','Vtheta Component (rad/s)','Velocity (mm/s)'] 

        iname = globalData[i]['name']
        print("Starting run " + iname + " ... this may take a while...")
        np.seterr(divide='ignore',invalid='ignore') # Divide 0 error protection

        # Converting to cylindrical coordinates.
        centreX = np.nanmean(posX)
        centreZ = np.nanmean(posZ)
        posX = posX - centreX
        posZ = posZ - centreZ
        posR = np.sqrt(posX**2 + posZ**2)
        posTheta = np.arctan2(posZ,posX)
        # NB posY is the same.

        # Setting up empty velocity arrays
        vR = np.full(len(posX), np.NaN)
        vY = np.full(len(posY), np.NaN)
        vTheta = np.full(len(posZ), np.NaN)

        # Linear regression - velocity calculations
        for i in range(jValue, len(time) - jValue):
            # Mask to avoid NaNs
            mask = ~np.isnan(posR[i - jValue:i + jValue]) & ~np.isnan(posY[i - jValue:i + jValue]) & ~np.isnan(posTheta[i - jValue:i + jValue])        

            if np.sum(mask) >= 2:
                # Defining Polynomials
                polyR = np.polyfit(time[i - jValue:i + jValue][mask], posR[i - jValue:i + jValue][mask], 1)
                polyY = np.polyfit(time[i - jValue:i + jValue][mask], posY[i - jValue:i + jValue][mask], 1)
                polyTheta = np.polyfit(time[i - jValue:i + jValue][mask], posTheta[i - jValue:i + jValue][mask], 1)

                # We want the gradient, which is the velocity.
                vR[i] = polyR[0]
                vY[i] = polyY[0]
                vTheta[i] = polyTheta[0]

        # Dividing by tip-speed, if necessary:
        if norm.lower() == 'tipspeed':
            vR = np.divide(vR,tipSpeed)
            vY = np.divide(vY,tipSpeed)
            vTheta = np.divide(vTheta, RPM*2*np.pi)

        # Trimming off outrageous velocities:
        vR = np.where(abs(vR) > vMax, np.nan, vR)
        vY = np.where(abs(vY) > vMax, np.nan, vY)
        vR = np.where(abs(posR) > TankD/2, np.nan, vR)
        vY = np.where(abs(posR) > TankD/2, np.nan, vY)
        posR = np.where(abs(posR) > TankD/2, np.nan, posR)
        posY = np.where(abs(posR) > TankD/2, np.nan, posY)
        vTheta = np.where(abs(vTheta) > 0.01, np.nan, vTheta)

        ### PLOTTING ###

        vMag = np.sqrt(vR**2 + vY**2 + vTheta**2)
        
        # Plot 1 - Distributions
        ax0[0].hist(vR, bins=bins, density=True, alpha=chartAlpha, label=iname)
        ax0[1].hist(vY, bins=bins, density=True, alpha=chartAlpha, label=iname)
        ax0[2].hist(vTheta, bins=bins, density=True, alpha=chartAlpha, label=iname)
        ax0[3].hist(vMag, bins=bins, density=True, alpha=chartAlpha, label=iname)

        # Plot Aesthetics
        for ax in ax0:
            ax.set_facecolor('white')  
        fig0.suptitle("Tank Cylindrical Velocity Distributions", fontsize = 14)

        ax0[0].set(xlabel=xlab[0], ylabel='Probability Density')
        ax0[1].set(xlabel=xlab[1], ylabel='Probability Density')
        ax0[2].set(xlabel=xlab[2], ylabel='Probability Density')
        ax0[3].set(xlabel=xlab[3], ylabel='Probability Density')

        ### PLOTTING PLOT 2 - AZIMUTHALLY AVERAGED VECTOR FIELDS ###
        print(np.nanmax(posR))
        print(np.nanmin(posR))
        print(np.nanmax(posY))
        print(np.nanmin(posY))
        
        # Axes labels
        xlab = '$(r/R)$'
        ylab = '$(h/H)$'

        # No cells in grid
        nr = int(nCellsWidth)
        cellSize = (np.nanmax(posR)/nr)*1.005 # Slightly enlarged to avoid issue with largest position points
        ny = int(np.ceil(np.nanmax(posY) / (cellSize)))

        # Creating grids
        vrGrid = np.zeros((ny, nr))
        vyGrid = np.zeros((ny, nr))
        nGrid = np.zeros((ny, nr))

        # For every particle position and velocity data point
        for i in range(len(vR)):
            r, y, vr, vy = posR[i], posY[i], vR[i], vY[i]

            # Exclude data points with NaN values
            if not np.isnan(r) and not np.isnan(y) and not np.isnan(vr) and not np.isnan(vy):
                ir = int(r / (cellSize))
                iy = int(y / (cellSize))

                nGrid[iy, ir] += 1  # Add one to the count
                vrGrid[iy, ir] += vr  # Adding the velocity to the velocity count
                vyGrid[iy, ir] += vy

        # To remove spurious velocities
        vrGrid[nGrid <= spuriousFilter] = np.nan
        vyGrid[nGrid <= spuriousFilter] = np.nan

        # Averaging Velocities
        vrGrid = np.nan_to_num(vrGrid / nGrid)
        vyGrid = np.nan_to_num(vyGrid / nGrid)

        # Defining position grids for plotting
        scaleR, scaleY = np.meshgrid(
            np.linspace(np.nanmin(posR)/(0.5*TankD), np.nanmax(posR)/(0.5*TankD), nr),  
            np.linspace(np.nanmin(posY)/(TankD), np.nanmax(posY)/(TankD), ny))

        # Normaise the plot vectors?
        if normVectors == True:
            normGrid = np.sqrt(vrGrid ** 2 + vyGrid ** 2)
        else:
            normGrid = 1

        # Define figure
        fig1, ax1 = plt.subplots(figsize=(6, 6), tight_layout=True)
        ax1.quiver(scaleR, scaleY, vrGrid / normGrid, vyGrid / normGrid, color='black', width=vWidth)

        # Plot Aesthetics
        rows, cols = nGrid.shape
        AR = 2

        ax1.set(ylabel=ylab, xlabel=xlab)
        ax1.set_aspect(AR)
        ax1.grid(False)

        ax1.set_xlabel(xlab, fontsize=17)
        ax1.set_ylabel(ylab, fontsize=17)
        ax1.tick_params(axis='x', labelsize=16)  # Set font size for x-axis tick labels
        ax1.tick_params(axis='y', labelsize=16)

        ax1.set_facecolor('white')
        if addTitle == True:
            fig1.suptitle("Velocity Field: " + iname, fontsize=14)
            
    ### OUTSIDE LOOP MISC ###
        
    # Set labels and legends for Plot 1
    for ax in ax0:
        ax.set_facecolor('white')  
        ax.legend()




