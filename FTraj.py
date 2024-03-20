## Filename: FTraj.py
## Date: 30/12/23
## Author: Itelcontar
## Description: Python3 module for creating basic plots of trajectories.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. length - Number of rows to use for trajectory plotting.
# 3. angle - Integer numpy array of same dimensions as exp. Values either 1,2 or 3. 1 is X-Y plane, 2 is Z-Y plane, and 3 is X-Z plane. Note PEPT covention holds Y as vertical, not Z.
    
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from FGlobalSettings import globalData

### FUNCTION ###

def fTraj(exp, length, angle):
    # Predefining figure
    fig, ax = plt.subplots(figsize=(6, 6))

    # For each experiment specified
    for i in range(0,len(exp)):
        # Collecting data
        posX = globalData[exp[i]]['posX']
        posY = globalData[exp[i]]['posY']
        posZ = globalData[exp[i]]['posZ']
        vX = globalData[exp[i]]['vX']
        name = globalData[exp[i]]['name']
        
        # Angle identification and adjustment
        if angle[i] == 2:
            posX = posZ
            xlab = 'z $(mm)$'
            ylab = 'y $(mm)$'
        elif angle[i] == 3:
            posY = posZ
            xlab = 'X $(mm)$'
            ylab = 'z $(mm)$'
        else:
            xlab = 'x $(mm)$'
            ylab = 'y $(mm)$'

        # Disallowing NaN rows
        NaNMask = np.full(len(posX),True)
        NaNMask = np.logical_not(np.isnan(vX)) 
        posMask = np.logical_not(np.isnan(posX))
        NaNMask = np.logical_and(posMask,NaNMask) 

        # Finding a random range
        sampleRange = []
        randomInt = 0
        n = 0
        
        # Use a while loop to generate a random range until it contains no NaN values
        while np.sum(np.isnan(posX[randomInt:randomInt + length])) > 0 or len(sampleRange) == 0:
            n = n + 1
            randomInt = np.random.randint(0, len(posX) - length + 1)
            sampleRange = posX[randomInt:randomInt + length]
            if n >1000:     # Breakout clause to prevent infinite loop
                length = round(length/2)
                n = 0
                if length == 1:
                    print("Length dropped to 0. Too many NaNs for plotting.")    
                    break
                  
        # Plotting using MatPlotLib
        linePos, = ax.plot(posX[randomInt:randomInt + length], posY[randomInt:randomInt + length], label=name)
        ax.set(xlabel=xlab, ylabel=ylab)
        ax.set_title("Particle Position in 2D Space: ", fontsize=14)
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend()
        
