## Filename: fVDistribSlice.py
## Date: 17/01/24
## Author: Itelcontar
## Description: Python3 module for creating Velocity Distrib Slice Plotting (i.e. distributions at different slices of height).
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. angle - Integer value. Values either 1,2 or 3. 1 is X-Y plane, 2 is Z-Y plane, and 3 is X-Z plane. Note PEPT covention holds Y as vertical, not Z.
# 3. segments - Number of slices to include in histogram plots.
# 4. velocity Spec - Which velocity component distribution to create? String: 'v', 'vX', 'vY', 'vZ'
# 5. bins - No. bins per distribution. Integer.
# 6. norm - Normalise velocities to tank diameter? Boolean True/False.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from FGlobalSettings import globalData
import matplotlib.ticker as mtick

### FUNCTION ###

def fVDistribSlice(exp, angle, segments, velocitySpec, bins, norm):

    # Initializing figure with GridSpec
    fig = plt.figure(figsize=(8, 2 * segments))
    gs = GridSpec(segments, 2, width_ratios=[0.03, 0.97], hspace=0.25)

    # Create the continuous axis on the left (plot part 2)
    leftAx = plt.subplot(gs[:, 0])
    leftAx.set_ylim(0, segments)

    # Set the number of ticks on the continuous left axis
    numTicks = segments + 1
    leftAx.set_yticks(np.linspace(0, segments, numTicks))

    # Initialize heightRange before the loop
    heightRange = []

    # Loop for every experiment in exp
    for n in exp:
        if norm == True:  
            # Selecting Velocities from Dictionary and normalizing. 
            posX = globalData[n]['posX'] / globalData[n]['TankD']
            posY = globalData[n]['posY'] / globalData[n]['TankD']
            posZ = globalData[n]['posZ'] / globalData[n]['TankD']
            vX = globalData[n]['vX'] / globalData[n]['TankD']
            vY = globalData[n]['vY'] / globalData[n]['TankD']
            vZ = globalData[n]['vZ'] / globalData[n]['TankD']
            title = "Tank Velocity (Normalised) Probability Density Distributions - Normalised with Tank Diameter"
        else:  
            # Selecting Velocities from Dictionary and normalizing only the positions. 
            posX = globalData[n]['posX'] / globalData[n]['TankD']
            posY = globalData[n]['posY'] / globalData[n]['TankD']
            posZ = globalData[n]['posZ'] / globalData[n]['TankD']
            vX = globalData[n]['vX']
            vY = globalData[n]['vY']
            vZ = globalData[n]['vZ']
            title = "Tank Velocity (Abs) Probability Density Distributions - Normalised with Tank Diameter"

        # Now we consider what velocity component we want for the plot
        if velocitySpec == "v":
            V = np.sqrt(vX**2 + vY**2 + vZ**2)
            if norm == True:
                xlab = "Velocity Magnitudes (tankD/s)"
            else:
                xlab = "Velocity Magnitudes (mm/s)"
        elif velocitySpec == "vX":
            V = vX
            if norm == True:
                xlab = "X Component Velocities (tankD/s)"
            else:
                xlab = xlab = "X Component Velocities (mm/s)"
        elif velocitySpec == "vY":
            V = vY
            if norm == True:
                xlab = "Y Component Velocities (tankD/s)"
            else:
                xlab = xlab = "Y Component Velocities (mm/s)"
        elif velocitySpec == "vZ":
            V = vZ
            if norm == True:
                xlab = "Z Component Velocities (tankD/s)"
            else:
                xlab = xlab = "Z Component Velocities (mm/s)"
            
        # We need to check and change angle:
        if angle == 3:
            posY = posZ  # Top-down view
            vY = vZ
            ylab = 'Z Position $(d/D)$'
        else:
            # Do nothing
            ylab = 'Y Position $(h/H)$'
        
        # Update heightRange based on the current experiment data
        heightRange = np.linspace(np.nanmin(posY), np.nanmax(posY), segments + 1)
        nonNanMask = ~np.isnan(V)
        
        # Subplots with shared x-axis
        for i in range(segments):
            subAx = plt.subplot(gs[i, 1])
            iname = globalData[n]['name']

            # Finding our relevant Velocities and Plotting
            sliceData = V[(posY >= heightRange[i]) & (posY < heightRange[i + 1]) & nonNanMask]
            subAx.hist(sliceData, bins=bins, density=True, alpha=1/len(exp), label =iname) #weights=100*np.ones(len(sliceData)) / len(sliceData)
            
            # Aesthetics
            subAx.yaxis.tick_right()
            subAx.grid(True, linestyle='--', alpha=0.7)
            subAx.set_facecolor('#f0f0f0')
            #subAx.yaxis.set_major_formatter(mtick.PercentFormatter(100))
            #subAx.set_ylim(0, 1)
            if i != segments-1:
                subAx.set_xticklabels([])

    # Plot Aesthetics for the last subplot
    subAx.set(xlabel=xlab)
    # And the first
    plt.subplot(gs[0, 1]).legend()
    plt.subplot(gs[0, 1]).set_title(title)

    # For Second Plot, which really acts as an axis.
    leftAx.set_ylabel(ylab)
    leftAx.set_xticks([]) # Constraining to make thin.
    leftAx.grid(True, linestyle='--', alpha=0.7)
    leftAx.set_facecolor('#f0f0f0')  # Light grey background
