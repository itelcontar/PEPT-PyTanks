## Filename: fVFieldDiff.py
## Date: 30/12/23
## Author: Itelcontar
## Description: Python3 module for creating Velocity Field Differential Plotting.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. noCellWidth - Integer number of cells for X-Axis. Other axes determiend by square pixel coniditon.
# 3. angle - Integer. Values either 1,2 or 3. 1 is X-Y plane, 2 is Z-Y plane, and 3 is X-Z plane. Note PEPT covention holds Y as vertical, not Z.
# 4. spuriousFilter - Integer, specifying minimum number of velocity data points in a pixel for the pixel to be included in the final plot. Removes 'spurious' velocities, with insufficient data for accurate averaging, or valid assumption of ergodicity.
# 5. plots - Boolean np array indicating which plots to create. [Vector Field, Magnitude error field, Contour Field, Angle error field] 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from FGlobalSettings import globalData

### FUNCTION ###

def fVFieldDiff(exp, noCellWidth, angle, spuriousFilter, plots):

    # First, collect some data
    time1 = globalData[exp[0]]['time']
    posX1 = globalData[exp[0]]['posX']
    posY1 = globalData[exp[0]]['posY']
    posZ1 = globalData[exp[0]]['posZ']
    vX1 = globalData[exp[0]]['vX']
    vY1 = globalData[exp[0]]['vY']
    vZ1 = globalData[exp[0]]['vZ']
    TankD1 = globalData[exp[0]]['TankD']
    time2 = globalData[exp[1]]['time']
    posX2 = globalData[exp[1]]['posX']
    posY2 = globalData[exp[1]]['posY']
    posZ2 = globalData[exp[1]]['posZ']
    vX2 = globalData[exp[1]]['vX']
    vY2 = globalData[exp[1]]['vY']
    vZ2 = globalData[exp[1]]['vZ']
    TankD2 = globalData[exp[1]]['TankD']
    
    name1 = globalData[exp[0]]['name']
    name2 = globalData[exp[1]]['name']
    
    # Next, normalise:
    vX1 = np.divide(vX1,TankD1)
    vY1 = np.divide(vY1,TankD1)
    vZ1 = np.divide(vZ1,TankD1)
    vX2 = np.divide(vX2,TankD2)
    vY2 = np.divide(vY2,TankD2)
    vZ2 = np.divide(vZ2,TankD2)
    posX1 = np.divide(posX1,0.5*TankD1)
    posY1 = np.divide(posY1,TankD1)
    posZ1 = np.divide(posZ1,0.5*TankD1)
    posX2 = np.divide(posX2,0.5*TankD2) 
    posY2 = np.divide(posY2,TankD2) 
    posZ2 = np.divide(posZ2,0.5*TankD2)

    # Centering:
    posX1 = posX1 - np.nanmean(posX1)
    posY1 = posY1 - np.nanmin(posY1)
    posZ1 = posZ1 - np.nanmean(posZ1)
    posX2 = posX2 - np.nanmean(posX2)
    posY2 = posY2 - np.nanmin(posY2)
    posZ2 = posZ2 - np.nanmean(posZ2)

    # Remove any datapoints outside the tank:
    posX1[posX1**2 + posZ1**2 > (0.505*TankD1)**2] = np.nan
    posZ1[posX1**2 + posZ1**2 > (0.505*TankD1)**2] = np.nan
    posX2[posZ2**2 + posX2**2 > (0.505*TankD2)**2] = np.nan
    posZ2[posZ2**2 + posX2**2 > (0.505*TankD2)**2] = np.nan
    
    # Check and change angle:
    if angle == 2:
        posX1 = posZ1  # Update variable names
        posX2 = posZ2
        vX1 = vZ1
        vX2 = vZ2
        xlab = '$(r/R)$'
        ylab = '$(h/H)$'
        yisy = True
    elif angle == 3:
        posY1 = posZ1  # Update variable names
        posY2 = posZ2
        vY1 = vZ1
        vY2 = vZ2
        xlab = '$(r/R)$'
        ylab = '$(r/R)$'
        yisy = False
    else:
        # Do nothing
        xlab = '$(r/R)$'
        ylab = '$(h/H)$'
        yisy = True

    # Then find a mutually beneficial grid shape
    nx = noCellWidth
    cellSize = [((np.nanmax(posX1) - np.nanmin(posX1))/nx)+0.0002, (((np.nanmax(posX2)) - (np.nanmin(posX2)))/nx) + 0.0002]
    
    # Aspect Ratio, to define square cells:
    if yisy == True:
        Aspect = (max(np.nanmax(posX1),np.nanmax(posX2)) - min(np.nanmin(posX1),np.nanmin(posX2)) / (max(np.nanmax(posY1),np.nanmax(posY2)) - min(np.nanmin(posY1),np.nanmin(posY2))))
    else:
        Aspect = 1
        
    # Finding number of cells for grids
    ny1 = int(np.ceil((np.nanmax(posY1) - np.nanmin(posY1)) / cellSize[0]) * Aspect) 
    ny2 = int(np.ceil((np.nanmax(posY2) - np.nanmin(posY2)) / cellSize[1]) * Aspect) 
    ny = max(ny1,ny2)
    
    def calculateVelocityField(time, posX, posY, posZ, vX, vY, vZ, cellSize, nx, ny, Aspect, yisy):

        # Creating Grids - All same
        count = 0
        vxgrid = np.zeros((ny, nx))
        vygrid = np.zeros((ny, nx))
        ngrid = np.zeros((ny, nx))
        xMin = np.nanmin(posX)
        if yisy == True:
            yMin = 0 # Should be 0?
        else:
            yMin = np.nanmin(posY)

        for i in range(len(vX)):
            x, y, vx, vy = posX[i], posY[i], vX[i], vY[i]

            # Exclude data points with NaN values
            if not np.isnan(x) and not np.isnan(y) and not np.isnan(vx) and not np.isnan(vy):
                ix = int((x - xMin)/ cellSize)
                iy = int(((y - yMin) / cellSize)*Aspect)

                # Ideally, we wouldn't need this if, but depending on how well your data is trimmed to the tank size, it may be useful to have
                if ix >= nx or iy >= ny: 
                    #Do nothing
                    count = count + 1
                else:
                    ngrid[iy, ix] += 1  # Add one to the count
                    vxgrid[iy, ix] += vx  # Adding the velocity to the velocity count
                    vygrid[iy, ix] += vy

        # Print how many velocities were out of range
        if count >10 :
            print("No. of Out of range Velocities: " + str(count))
            print("If >100: Check for particle time outside tank volume, and trim data.")

        # To remove spurious velocities
        vxgrid[ngrid <= spuriousFilter] = np.nan
        vygrid[ngrid <= spuriousFilter] = np.nan

        # Averaging velocities
        vxgrid = np.nan_to_num(vxgrid / ngrid)  
        vygrid = np.nan_to_num(vygrid / ngrid)

        norm_on = False 
        return vxgrid, vygrid
        # END of General Function

    # Calculate velocity fields for both sets of data
    vxgrid1, vygrid1 = calculateVelocityField(time1, posX1, posY1, posZ1, vX1, vY1, vZ1, cellSize[0], nx, ny, Aspect,yisy)
    vxgrid2, vygrid2 = calculateVelocityField(time2, posX2, posY2, posZ2, vX2, vY2, vZ2, cellSize[1], nx, ny, Aspect,yisy)

    if yisy == True:
        ymin = min(np.nanmin(posY1), np.nanmin(posY2))
        ymax = 1
        ydiff = (ymax-ymin)/ny
        AR = 2 # Aspect Ratio
    else:
        ymin = ymin = min(np.nanmin(posY1), np.nanmin(posY2))
        ymax = 1
        ydiff = (ymax-ymin)/ny
        AR = 1 # Aspect Ratio

    # X Limits are always the same, as X is always the same
    xMin = min(np.nanmin(posX1), np.nanmin(posX2))
    xMax = max(np.nanmax(posX1), np.nanmax(posX2))
    xDiff = (max(np.nanmax(posX1), np.nanmax(posX2))-min(np.nanmin(posX1), np.nanmin(posX2)))/nx
        
    # Making grids
    scale_x, scale_y = np.meshgrid(
        np.arange(xMin, xMax, xDiff), 
        np.arange(ymin, ymax, ydiff) 
        )

    ### PLOTTING ###
    
    if plots[1] == True:
        # Calculate the magnitude of vectors for both sets of data
        magnitude1 = np.sqrt(vxgrid1**2 + vygrid1**2)
        magnitude2 = np.sqrt(vxgrid2**2 + vygrid2**2)

        # Calculate the difference in magnitude
        magnitude_difference = magnitude1 - magnitude2

        # Plot the magnitude difference with imshow
        fig, ax = plt.subplots(figsize=(8, 8), tight_layout=True)

        # Plot the color-coded representation of the magnitude difference (imshow)
        im = ax.imshow(magnitude_difference, extent=[xMin,xMax,ymin,ymax], cmap='jet')

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Vector Magnitude Difference')

        # Plot Aesthetics
        ax.set(ylabel=ylab, xlabel=xlab)
        ax.set_aspect(1)
        ax.grid(True, linestyle='--', alpha=0.7)

        ax.set_facecolor('white')
        fig.suptitle("Magnitude Difference between Velocity Fields: " + name1 + " - " + name2, fontsize=14)
        ax.set_aspect(AR)
        
    if plots[0] == True:
        # Plot individual velocity fields on the same plot
        fig1, ax1 = plt.subplots(figsize=(7, 7), tight_layout=True)

        # Plot for vxgrid1 and vygrid1 in blue
        ax1.quiver(scale_x, scale_y, vxgrid1, vygrid1, color='black', width=0.003, label=name1)

        # Plot for vxgrid2 and vygrid2 in green
        ax1.quiver(scale_x, scale_y, vxgrid2, vygrid2, color='b', width=0.003, label=name2)

        ax1.set(ylabel=ylab, xlabel=xlab)
        ax1.set_aspect(AR)
        ax1.set_facecolor('white')

        ax1.set_xlabel(xlab, fontsize=17)
        ax1.set_ylabel(ylab, fontsize=17)
        ax1.tick_params(axis='x', labelsize=16)  # Set font size for x-axis tick labels
        ax1.tick_params(axis='y', labelsize=16)
        legend = ax1.legend(fontsize=16,loc='upper right')
        legend.get_frame().set_edgecolor('black')
        #fig1.suptitle(name1 + " & " + name2, fontsize=16)
        
    if plots[2] == True:
        # Defining Figure
        fig, ax = plt.subplots(figsize=(8, 8), tight_layout=True)

        # Stream plotting
        stream1 = ax.streamplot(scale_x, scale_y, vxgrid1, vygrid1, color='black')
        stream2 = ax.streamplot(scale_x, scale_y, vxgrid2, vygrid2, color='b')

        # Create proxy artists for the legend
        legend_labels = [name1, name2]
        legend_elements = [Patch(facecolor='black', label=legend_labels[0]),
                           Patch(facecolor='b', label=legend_labels[1])]

        # Plot Aesthetics
        ax.set(ylabel=ylab, xlabel=xlab)
        ax.set_aspect(AR)

        ax.set_facecolor('white')
        fig.suptitle("Velocity Field Contours: " + name1 + " & " + name2, fontsize=14)
        ax.set_xlabel(xlab, fontsize=17)
        ax.set_ylabel(ylab, fontsize=17)
        ax.tick_params(axis='x', labelsize=16)  # Set font size for x-axis tick labels
        ax.tick_params(axis='y', labelsize=16)
        ax.legend(handles=legend_elements, fontsize=16, loc='upper right')

        # Add a bit of whitespace around the plot
        ax.set_xlim(xMin*1.1, xMax*1.1)
        ax.set_ylim(ymin*1.1, ymax*1.1)
        
    if plots[3] == True:
        # Calculate the angle of vectors for both sets of data
        angle1 = np.arctan2(vygrid1, vxgrid1)
        angle2 = np.arctan2(vygrid2, vxgrid2)

        # Calculate the difference in angles
        angleDifference = np.degrees(np.abs(angle1 - angle2))

        # Plot the magnitude difference with imshow
        fig, ax = plt.subplots(figsize=(8, 8), tight_layout=True)

        # Plot the color-coded representation of the magnitude difference (imshow)
        im = ax.imshow(angleDifference, extent=[-1,1,ymin,ymax], cmap='jet')

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Vector Angle Difference')

        # Plot Aesthetics
        ax.set(ylabel=ylab, xlabel=xlab)
        ax.set_aspect(1)
        ax.grid(True, linestyle='--', alpha=0.7)

        fig.suptitle("Angle Difference between Velocity Fields: " + name1 + " - " + name2, fontsize=14)
        ax.set_aspect(AR)

        
