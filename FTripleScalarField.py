## Filename: FTripleScalarField.py
## Date: 22/02/24
## Author: Itelcontar
## Description: Python3 module for creating Velocity Scalar Fields, Residency Scalar Fields, or Velocity Quiver Plots, in a pretty format.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. cellSize - In mm, single dimension size for cubic voxels.
# 3. spuriousFilter - Integer, specifying minimum number of velocity data points in a voxel for the voxel to be included in the final plot. Removes 'spurious' velocities, with insufficient data for accurate averaging, or valid assumption of ergodicity.
# 4. plots - Boolean numpy array for generating various plots. [Velocity Scalar, Residency Scalar, Velocity Vector & Scalar].

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData

### FUNCTION ###

def fTripleScalarField(exp, cellSize, spuriousFilter, plots):
    # For speed, defining internal function for plotting individual scalar fields
    def scalarVPlot(exp, cellSize, angle, spuriousFilter, fig, ax):
        # Gathering data, NB all data is normalised with TankD
        TankD = globalData[exp]['TankD']
        posX = np.divide(globalData[exp]['posX'],TankD)
        posY = np.divide(globalData[exp]['posY'],TankD)
        posZ = np.divide(globalData[exp]['posZ'],TankD)
        vX = globalData[exp]['vX']
        vY = globalData[exp]['vY']
        vZ = globalData[exp]['vZ']
        name = globalData[exp]['name']
        vMag = np.sqrt(vX**2 + vY**2 + vZ**2)
        
        # Check and change angle
        if angle == 2:
            posX = posZ  # Alternative side view
            vX = vZ
            xlab = '$(r/R)$'
            ylab = '$(h/H)$'
            yisy = True
        elif angle == 3:
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

        # Defining no. cells in Grid
        nx = int(np.ceil(np.nanmax(posX)*TankD / (cellSize)))
        ny = int(np.ceil(np.nanmax(posY)*TankD / (cellSize)))

        # Creating grids, for velocities and number of the same
        scalGrid = np.zeros((ny, nx))
        ngrid = np.zeros((ny, nx))

        # Running through every velocity data point
        for i in range(len(vX)):
            x, y, v = posX[i], posY[i], vMag[i]

            # Exclude data points with NaN values
            if not np.isnan(x) and not np.isnan(y) and not np.isnan(v):
                ix = int(x *TankD/ (cellSize))
                iy = int(y *TankD/ (cellSize))

                # Having ensured value is usable, it is recorded by adding to the grids
                ngrid[iy, ix] += 1  
                scalGrid[iy, ix] += v  # Adding the velocity to the velocity grid           

        # To remove spurious velocities
        scalGrid[ngrid <= spuriousFilter] = np.nan

        # Averaging, to find real velocities, not sum
        scalGrid = np.nan_to_num(scalGrid / ngrid)

        # Before calculating position grids, understand axes limits, and centre on 0
        # Before calculating position grids, need to know plot type (xy/zy or z/x) and set appropriate limits
        if yisy == True:
            ymin = 0
            ymax = 1
            ydiff = 1/ny
            AR = 2 # Aspect Ratio
        else:
            ymin = -1
            ymax = 1
            ydiff = 2/ny
            AR = 1 # Aspect Ratio

        # Plotting
        im = ax.imshow(scalGrid, origin = 'lower', extent=[-np.nanmax(posX),np.nanmax(posX),ymin,ymax], cmap = 'plasma')#, vmin=0, vmax=650)

        # Plot Aesthetics
        ax.set_ylabel(ylab, fontsize = 11)
        ax.set_xlabel(xlab, fontsize = 11)
        #ax.tick_params(axis='x', labelsize=14)  # Set font size for x-axis tick labels
        #ax.tick_params(axis='y', labelsize=14)
        ax.set_aspect(AR)
        return im 


### ### ### ### ### ### ### ### ### #### ### ### #### ### ## ### ### ## ## #### ### ## ## ### ## ### #


    def scalarTPlot(exp, cellSize, angle, fig, ax):
        # Data gathering 
        TankD = globalData[exp]['TankD']
        posX = np.divide(globalData[exp]['posX'],TankD)
        posY = np.divide(globalData[exp]['posY'],TankD)
        posZ = np.divide(globalData[exp]['posZ'],TankD)
        time = globalData[exp]['time']
        name = globalData[exp]['name']
        cellSize = cellSize / TankD
        
        # First check and change angle
        if angle == 2:
            posX = posZ  # Alternative side view
            xlab = '$(r/R)$'
            ylab = '$(h/H)$'
            yisy = True
        elif angle == 3:
            posY = posZ  # Top down view
            xlab = '$(r/R)$'
            ylab = '$(r/R)$'
            yisy = False
        else:
            # Do nothing
            xlab = '$(r/R)$'
            ylab = '$(h/H)$'
            yisy = True

        posMid = np.zeros((len(posX) - 1, 3))
        dt = np.zeros(len(posX) - 1)

        for i in range(0,len(posX)-1): # Finding average positions (effectively cleans data a tad)
            posMid[i,0] = (posX[i] + posX[i+1]) / 2 
            posMid[i,1] = (posY[i] + posY[i+1]) / 2
            posMid[i,2] = (posZ[i] + posZ[i+1]) / 2
            dt[i] = time[i+1] - time[i]

        # No cells in grid
        nx = int(np.ceil(np.nanmax(posX) / cellSize))
        ny = int(np.ceil(np.nanmax(posY) / cellSize))
        
        # Creating grid
        Resgrid = np.zeros((ny,nx))

        # Mask check to find and ignore NaNs
        ResMask = np.logical_not(np.logical_or( np.isnan(posX), np.isnan(posY)))

        # Loop through and add up time spent in each cell
        for i in range(0,len(posX) - 1):
            if ResMask[i] == True: # NB true means we have an clean data point with no NaNs
                ix = int(posX[i] / cellSize)
                iy = int(posY[i] / cellSize)
                Resgrid[iy,ix] = Resgrid[iy,ix] + dt[i]  

        norm_on = True # Handle

        if norm_on == True:
            norm = Resgrid / np.sum(Resgrid)
        else:
            norm = 1
        
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

        # Plotting
        im = ax.imshow(Resgrid, origin = 'lower', extent=[-np.nanmax(posX),np.nanmax(posX),ymin,ymax], cmap = 'plasma')#, vmin=0, vmax=0.5)

        # Plot Aesthetics
        ax.set(ylabel=ylab, xlabel=xlab)

        rows, cols = Resgrid.shape
        AR = cols/rows
        ax.set_aspect(AR)
        return im
    
### ### ### ### ### ### ### ### ### #### ### ### #### ### ## ### ### ## ## #### ### ## ## ### ## ### #

    def vectorVPlot(exp, cellSize, angle, spuriousFilter, fig, ax):
        # Data gathering, NB Normalised to TankD by default
        TankD = globalData[exp]['TankD']
        posX = np.divide(globalData[exp]['posX'],TankD)
        posY = np.divide(globalData[exp]['posY'],TankD)
        posZ = np.divide(globalData[exp]['posZ'],TankD)
        vX = globalData[exp]['vX']
        vY = globalData[exp]['vY']
        vZ = globalData[exp]['vZ']
        name = globalData[exp]['name']

        # First check and change angle
        if angle == 2:
            posX = posZ  # Alternative side view
            vX = vZ
            xlab = '$(r/R)$'
            ylab = '$(h/H)$'
            yisy = True
        elif angle == 3:
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

        # Creating grids for velocity components, and a number count grid
        vxgrid = np.zeros((ny, nx))
        vygrid = np.zeros((ny, nx))
        ngrid = np.zeros((ny, nx))

        # For each velocity data point
        for i in range(len(vX)):
            x, y, vx, vy = posX[i], posY[i], vX[i], vY[i]

            # Exclude data points with NaN values
            if not np.isnan(x) and not np.isnan(y) and not np.isnan(vx) and not np.isnan(vy):
                ix = int(x *TankD/ (cellSize))
                iy = int(y *TankD/ (cellSize))

                ngrid[iy, ix] += 1  # Add one to the count
                vxgrid[iy, ix] += vx  # Adding the velocity to the velocity count
                vygrid[iy, ix] += vy

        # To remove spurious velocities
        vxgrid[ngrid <= spuriousFilter] = np.nan
        vygrid[ngrid <= spuriousFilter] = np.nan

        # Averaging
        vxgrid = np.nan_to_num(vxgrid / ngrid)
        vygrid = np.nan_to_num(vygrid / ngrid)

        # Before calculating position grids, need to know plot type (xy/zy or z/x) and set appropriate limits
        if yisy == True:
            ymin = 0
            ymax = 1
            ydiff = 1/ny
            AR = 2 # Aspect Ratio
        else:
            ymin = -1
            ymax = 1
            ydiff = 2/ny
            AR = 1 # Aspect Ratio
            
        # Plotting
        scale_x, scale_y = np.meshgrid(
            np.arange(-1, 1, 2/nx),
            np.arange(ymin, ymax, ydiff)
        )
        
        # Plotting
        im = ax.quiver(scale_x, scale_y, vxgrid, vygrid, color = 'white',width = 0.003)
        ax.set_xlim(-1,1)
        ax.set_ylim(ymin,ymax)

### ### ### ### ### ### ### ### ### #### ### ### #### ### ## ### ### ## ## #### ### ## ## ### ## ### #

    # Now the functions are defined, a for loop is used to run through each experiment, generate all three plots, and put them together
            
    # Scalar plot for velocities
    if plots[0] == True:
        for n in exp:
            
            # Defining Big Figure
            fig, ax =  plt.subplots(1,3,figsize=(18, 5))
            name = globalData[n]['name']
            
            scalarVPlot(n, cellSize, 1, spuriousFilter, fig, ax[0])
            scalarVPlot(n, cellSize, 2, spuriousFilter, fig, ax[1])
            im = scalarVPlot(n, cellSize, 3, spuriousFilter, fig, ax[2])
            
            # Add a single colorbar for all subplots to the right
            cbar = fig.colorbar(im, ax=ax, orientation='vertical', fraction=0.04, pad=0.025)
            cbar.set_label("mm/s", fontsize=12)

            #Aesthetics:
            fig.suptitle(name + " Velocity Scalar Field", fontsize = 14)

    # Scalar plot for residence time
    if plots[1] == True:
        for n in exp:
            
            # Defining big figure
            fig, ax =  plt.subplots(1,3,figsize=(18, 5)) 
            name = globalData[n]['name']
            
            scalarTPlot(n, cellSize, 1, fig, ax[0])
            scalarTPlot(n, cellSize, 2, fig, ax[1])
            im = scalarTPlot(n, cellSize, 3, fig, ax[2])
            
            # Add a single colorbar for all subplots to the right
            cbar = fig.colorbar(im, ax=ax, orientation='vertical', fraction=0.04, pad=0.025)
            cbar.set_label("%t", fontsize=12)

            # Aesthetics
            fig.suptitle(name + " Residency Scalar Field", fontsize = 14)

    # Quiver plot for velocities, overlayed on scalar field
    if plots[2] == True: # figsize default 18,5
        for n in exp:

            # Defining big figure
            fig, ax =  plt.subplots(1,3,figsize=(18, 5))
            name = globalData[n]['name']
            
            scalarVPlot(n, cellSize, 1, spuriousFilter, fig, ax[0])
            scalarVPlot(n, cellSize, 2, spuriousFilter, fig, ax[1])
            im = scalarVPlot(n, cellSize, 3, spuriousFilter, fig, ax[2])
            
            # Add a single colorbar for all subplots to the right
            cbar = fig.colorbar(im, ax=ax, orientation='vertical', fraction=0.04, pad=0.025)
            cbar.set_label("mm/s", fontsize=12)

            # Now add a vector field on top
            vectorVPlot(n, cellSize*3, 1, spuriousFilter, fig, ax[0])
            vectorVPlot(n, cellSize*3, 2, spuriousFilter, fig, ax[1])
            vectorVPlot(n, cellSize*3, 3, spuriousFilter, fig, ax[2])

            # Aesthetics
            fig.suptitle(name + " Velocity Scalar & Vector Fields", fontsize = 14)
        

