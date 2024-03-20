## Filename: FVRad.py
## Date: 05/03/24
## Author: Itelcontar
## Description: Python3 module for plotting radial velocity distributions.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. component - Either 'v' or 'vy, indicates velocity component to study.

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData

### FUNCTION ###

def fVRad(exp,component):
    # Defining Figure
    fig, ax = plt.subplots(figsize=(6, 6), tight_layout=True)
    lineStyles = np.array(['None','-', '--', '-.', ':', ' ', ''])
    count = 0
    
    # Running through each experiment specified
    for i in range(0,len(exp)):
        # Data gathering
        TankD = globalData[exp[i]]['TankD']
        time = globalData[exp[i]]['time']
        posX = np.divide(globalData[exp[i]]['posX'],TankD)
        posY = np.divide(globalData[exp[i]]['posY'],TankD)
        posZ = np.divide(globalData[exp[i]]['posZ'],TankD)

        # Finding impeller tip speed
        radius = TankD / 2
        RPM = globalData[exp[i]]['RPM']
        tipSpeed = (1/3) * radius * 2 * np.pi * RPM / 60 # (1/3) * 2piR * RPM / 60 => mm/s tip speed            

        # Velocities
        vX = np.divide(globalData[exp[i]]['vX'],tipSpeed)
        vY = np.divide(globalData[exp[i]]['vY'],tipSpeed)
        vZ = np.divide(globalData[exp[i]]['vZ'],tipSpeed)
        name = globalData[exp[i]]['name']
        vMag = np.sqrt(vX**2 + vY**2 + vZ**2)
            
        # Axes labels
        xlab = '$(r/R)$'
        count = count + 1

        # Component
        if component.lower() == 'v':
            V = vMag
            ylab = '$(V/V_{tip})$'
        elif component.lower() == 'vy' or component.lower() == 'y':
            V = vY
            ylab = '$(V_{Y}/V_{tip})$'

        # Radial Position and grid
            # Normalise Position
        posX = 2* (posX - np.nanmean(posX)) # 2* from radius not diameter
        posZ = 2* (posZ - np.nanmean(posZ))
        posR = np.sqrt(posX**2 + posZ**2)
        nr = 49 # Change as desired
        cellSize = np.nanmax(posR)/nr        
        
        # Creating grids
        vGrid = np.zeros(nr+1)
        nGrid = np.zeros(nr+1)

        # Cycling through all velocities, and sorting into bins
        for i in range(len(vX)):
            r, v = posR[i], V[i]

            # Exclude data points with NaN values
            if not np.isnan(r) and not np.isnan(v):
                ir = int(r/(cellSize))

                nGrid[ir] += 1  # Add one to the count
                vGrid[ir] += v  # Adding the velocity to the velocity count

        vGrid = np.nan_to_num(vGrid / nGrid)  # Averaging
        R = np.linspace(0,np.nanmax(posR),nr+1)
        
        # Define figure
        ax.plot(R,vGrid,label = name, color = 'black', linestyle = lineStyles[count])

        ax.set(ylabel=ylab, xlabel=xlab)
        ax.grid(False)
        ax.set_xlim([0,1])
        ax.legend()
        ax.set_facecolor('white')
        #fig.suptitle("Radial Velocity Distribution: " + name, fontsize=14) 
