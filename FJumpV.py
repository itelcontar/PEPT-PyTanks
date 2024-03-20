## Filename: FJumpV.py
## Date: 03/02/24
## Author: Itelcontar
## Description: Python3 module for calculating jump velocities.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. mjt - Minimum jump time, time required to be classified as a jump (s).
# 3. mjd - Minimum jump distance, distance required to be classified as a jump (TankD).
# 4. plots - Boolean array indicated which plots to create [Jump Velocities, Jump Times, Jump Distances].

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData

### FUNCTION ###

def fJumpV(exp, mjt, mjd, plots):
    # Predefine plots
    fig, ax = plt.subplots(figsize=(8, 8), tight_layout=True)
    fig1, ax1 = plt.subplots(figsize=(8, 8), tight_layout=True)
    fig2, ax2 = plt.subplots(figsize=(8, 8), tight_layout=True)

    # Running through each specified experiment
    for n in range(0, len(exp)):
        TankD = globalData[exp[n]]['TankD']
        time = globalData[exp[n]]['time']
        posY = np.divide(globalData[exp[n]]['posY'], TankD)
        name = globalData[exp[n]]['name']

        # Initialization of Lists
        jumpTime = []
        jumpDistance = []
        jumpVelocities = np.array([], dtype=float)

        # Boolean variable for jumpSwitch
        jumpSwitch = False

        # Now calculate jump velocities by cycling through full trajectory
        for i in range(1, len(time) - 1):
            itime = time[i]
            posOld = posY[i - 1]
            posNew = posY[i]

            if posNew > posOld:  # Paricle moves upwards
                if not jumpSwitch:
                    jumpSwitch = True
                    jumpStartPos = posNew
                    jumpStartTime = itime
            else:
                if jumpSwitch:
                    jumpSwitch = False
                    if posNew - jumpStartPos >= mjd and itime - jumpStartTime >= mjt:
                        # A good jump was detected
                        jumpTime.append(itime - jumpStartTime)
                        jumpDistance.append(posNew - jumpStartPos)
                        jumpVelocities = np.append(jumpVelocities, (posNew - jumpStartPos) / (itime - jumpStartTime))

        if plots[0] == True:
            # Plotting
            ax.hist(jumpVelocities, bins=50, density=True, label=name, alpha=1 / len(exp))

            # Plot Aesthetics
            ax.set(ylabel="Probability", xlabel="Velocities (TankD/s)")
            ax.grid(False)
            fig.suptitle("Bubble Jump Velocities", fontsize=14)
            ax.legend()
            
        if plots[1] == True:
            # Plotting
            ax1.hist(jumpTime, bins=50, density=True, label=name, alpha=1 / len(exp))

            # Plot Aesthetics
            ax1.set(ylabel="Probability", xlabel="Jump Times (s)")
            ax1.grid(False)
            fig1.suptitle("Bubble Jump Times", fontsize=14)
            ax1.legend()
            
        if plots[2] == True:
            # Plotting
            ax2.hist(jumpDistance, bins=50, density=True, label=name, alpha=1 / len(exp))

            # Plot Aesthetics
            ax2.set(ylabel="Probability", xlabel="Jump Distances (TankD)")
            ax2.grid(False)
            fig2.suptitle("Bubble Jump Distances", fontsize=14)
            ax2.legend()
