## Filename: FRCTime.py
## Date: 31/01/24
## Author: Itelcontar
## Description: Python3 module for creating recirculation time Distributions.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. cellSize - Integer width no. of velocity field cells.
# 3. plane - Indicates the central plane of interest (typically impeller plane) in dimensionless height (i.e. 0.33).
# 4. plots - Boolean array indicated which plots to create [Recirculation time distributions, Scalar field for time taken to return to impeller plane, Example Loop plots]
# 5. mlt - Maximum loop time, in seconds. All loops larger than this value are excluded. 
# 6. noLoops - Number of loops to plot on example loop plot (plots[2])
# 7. PlaneTolerance - Dimensionless height baove and below plane of interest that particle must cross, to consider a loop 'good'.

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData

### FUNCTION ###

def fRCTime(exp, cellSize, plane, plots, mlt, noLoops,PlaneTolerance):
    fig1, ax1 = plt.subplots(figsize=(6, 6), tight_layout=True)

    # Initialising figure variables
    line_styles = ['-', '--', ':', '-', '--', '-.', ':']
    colours = ['black', 'black', 'black','red','red','red','red']
    count = -1

    # Running through each experiment specified
    for n in range(0, len(exp)):
        TankD = globalData[exp[n]]['TankD']
        time1 = globalData[exp[n]]['time']
        posX = np.divide(globalData[exp[n]]['posX'], 0.5*TankD)
        posY = np.divide(globalData[exp[n]]['posY'], TankD)
        posZ = np.divide(globalData[exp[n]]['posZ'], 0.5*TankD)
        name = globalData[exp[n]]['name']

        # Centering plot:
        posX = posX - np.nanmean(posX)
        posZ = posZ - np.nanmean(posZ)

        # Where is the impeller plane?
        pImp = plane  # Because we normalized positions. Note data must be trimmed to the bottom of tank @ y = 0.
        posYRel = posY - pImp

        # Define two more planes - symetric above and below the impeller plane
        checkPlaneAbove = PlaneTolerance
        checkPlaneBelow = -PlaneTolerance

        # 0. Finding first location above top plane tolerance 
        i = 1
        E1 = 0
        E2 = 0
        while E1*E2 >= 0: # Loop breaks at i = k+1 where k is index of first sign change 
            E2 = posYRel[i] - checkPlaneAbove
            E1 = posYRel[i-1] - checkPlaneAbove
            if E1 == np.nan:
                E1 = 0
            if E2 == np.nan:
                E2 = 0
            i = i + 1 # Add to the count
        k = i - 2 # Index of last point in first half-loop.
        
        # 1. Initialise new variables
        t2LoopEnd = np.zeros(len(posYRel))  # Time until end of loop for which a position is in.
        loopTimes = np.array([])            # List of total loop times.
        loopIndices = np.empty((0, 2))      # List of loop indices.           
        iLoopStart = k+1                    # Index of latest loop first point, starts with first full loop.
        writeDataSwitch = False
        badLoopSwitch = False
        loopType = ''
        crossings = 0                       # Count of no. crossings
        
        t2LoopEnd[0:k] = np.nan # Blanking first half loop out

        # 2. Loop through full set of position data from first loop:
        for i in range(k+1,len(posYRel)):
            # 2.1 Check if last loop or not:
            if i < len(posYRel) - 1: 
                # No - collecting data
                Y1 = posYRel[i] # Current Y
                Y2 = posYRel[i+1] # Next Y
            else:
                # Yes - last loop in data
                badLoopSwitch = True
                writeDataSwitch = True
            
            # 2.2 Loop Checking for NaNs
            if np.isnan(Y1*Y2) == np.nan: 
                badLoopSwitch = True

            # 2.3 Sign change over top plane? (I.e. has loop finished?)
            if (Y1-checkPlaneAbove)*(Y2-checkPlaneAbove) < 0:
                if crossings == 1:
                    writeDataSwitch = True
                    
                    # 2.4 Check if the loop has passed around the verification planes
                    loopPosY = posYRel[iLoopStart:i] # This is the list of Y positions in the loop.
                    if len(loopPosY) > 0:
                        if np.nanmax(loopPosY) > checkPlaneAbove and np.nanmin(loopPosY) < checkPlaneBelow:
                            # We have a good loop
                            loopType = 'good'
                            #print("Good loop")
                        else:
                            badLoopSwitch = True
                            #print("Bad loop")
                    else:
                        badLoopSwitch = True
                else:
                    crossings = crossings + 1

            # 2.4 Check if loop time is smaller than maximum loop time (mlt)
            if time1[i] - time1[iLoopStart] > mlt: 
                badLoopSwitch = True
            
            # 2.5. Data Writing and loop restart
            if writeDataSwitch == True: # If write data was flagged:
                if badLoopSwitch == True:
                    t2LoopEnd[iLoopStart:i] = np.nan # Wipe out the whole loop.
                    iLoopStart = i + 1 # New start index
                    badLoopSwitch = False # Reset the switches
                    writeDataSwitch = False
                    crossings = 0
                else:
                    t2LoopEnd[iLoopStart:i] = time1[i] - time1[iLoopStart:i] # Write in full loop of data
                    loopTime = time1[i] - time1[iLoopStart] # Calculate loop time
                    loopTimes = np.append(loopTimes, loopTime) # Add to archived loop times
                    loopIndices = np.vstack([loopIndices, np.array([iLoopStart, i])])
                    iLoopStart = i + 1 # New start index
                    writeDataSwitch = False # Reset switch
                    crossings = 0
        
        ### PLOT 1 ###   ---    Map of time to return to impeller plane
                    
        # No cells in Grid.
        nx = int(np.ceil(np.nanmax(posX) / cellSize))
        ny = int(np.ceil(np.nanmax(posY) / cellSize))

        # Creating Grid
        rcGrid = np.zeros((ny, nx))
        nGrid = np.zeros((ny, nx))

        # Mask Check, so we don't use NaNs!
        Mask = np.logical_not(np.logical_or(np.isnan(t2LoopEnd),np.isnan(posX)))

        if plots[0] == True:
            # Loop through and add up time until we reach the impeller plane.
            for i in range(0, len(time1) - 1):
                if Mask[i] == True:  # NB true means we have a clean data point with no NaNs.
                    ix = int(posX[i] / cellSize)
                    iy = int(posY[i] / cellSize)

                    # We find the time until we return to the impeller plane.
                    rc = t2LoopEnd[i]

                    # Add to grid.
                    rcGrid[iy, ix] = rcGrid[iy, ix] + rc
                    nGrid[iy, ix] = nGrid[iy, ix] + 1
            
            norm_on = True  # HANDLE
            if norm_on == True:
                nGrid = np.where(nGrid == 0, np.nan, nGrid)
                rcGrid = np.divide(rcGrid, nGrid)
            #rcGrid = np.where(rcGrid > 10, 10, rcGrid) # Removing excessively high values where particle is stuck

            # Plotting
            fig0, ax0 = plt.subplots(figsize=(6, 6), tight_layout=True)
            im = ax0.imshow(rcGrid, origin='lower', extent=[np.nanmin(posX), np.nanmax(posX), np.nanmin(posY), np.nanmax(posY)], cmap='jet')

            # Plot Aesthetics
            ax0.set(ylabel="h/H", xlabel="r/R")
            ax0.set_aspect(2)
            ax0.grid(False)
            fig0.suptitle("Time to Return to Impeller-Plane: " + name, fontsize=14)
            colorbar = plt.colorbar(im, ax=ax0)


        ### PLOT 2 ###   ---   Distribution of recirculation times.

        # Distribution Plot
        if plots[1] == True:
            n, bins, _ = ax1.hist(loopTimes, bins=100, density=True, alpha=0)
            count = count + 1

            # Calculate bin centers
            bin_centers = (bins[:-1] + bins[1:]) / 2
            
            # Calculate the weighted average of the probability density function
            weightedAverage = np.sum(bin_centers * n) / np.sum(n)

            # Print the weighted average
            print("Weighted Average:", weightedAverage)

            # Calculate cumulative distribution
            cumulative = np.cumsum(n)
            cumulative = np.divide(cumulative, np.nanmax(cumulative))

            # Plot the line using bin centers and histogram values
            ax1.plot(bin_centers, cumulative, label=name, linestyle=line_styles[count], color = colours[count])

            ax1.grid(True, linestyle='--', alpha=0.7)
            ax1.set_facecolor('white')
            ax1.legend()
            ax1.set_xlabel('Recirculation Time (s)', fontsize = 15)
            ax1.set_ylabel('Cumulative Probability', fontsize = 15)
            ax1.tick_params(axis='x', labelsize=13.5)  # Set font size for x-axis tick labels
            ax1.tick_params(axis='y', labelsize=13.5)
            fig1.suptitle("Recirculation Times - 2 Plane, Separation: " + str(2*PlaneTolerance))

        ### PLOT 3 ###   ---   Trajectories.
        
        # Example Loop Trajectories
        if plots[2] == True:
            # Create a New Plot for Each EXP.
            fig3, ax3 = plt.subplots(figsize=(6, 6), tight_layout=True)

            # Finding some indices for random loops
            indices = loopIndices[np.random.choice(loopIndices.shape[0], size=noLoops, replace=False)]
            
            # Now, similar to FTraj.py, plot the trajectories between the indices identified, to visualise the loops
            for i in range(1,noLoops):
                # Plot Loops
                iStart = int(indices[i-1,0])
                iEnd   = int(indices[i-1,1])
                sampleX = posX[iStart:iEnd]
                sampleY = posY[iStart:iEnd]
                linePos, = ax3.plot(sampleX, sampleY, color = 'blue', label='Recirculation Loops')

            # Plane Lines
            ax3.axhline(y=pImp, color = 'red', linestyle='--', label = 'Impeller Plane')
            ax3.axhline(y=(pImp + PlaneTolerance), color = 'black', linestyle='--', label = 'Verification Planes')
            ax3.axhline(y=(pImp - PlaneTolerance), color = 'black', linestyle='--')
            
            # Finishing touches
            ax3.set_title("2D Recirculations: " + name, fontsize=14)
            ax3.set_xlabel("r/R",fontsize = 14)
            ax3.set_ylabel("h/H",fontsize = 14)
            ax3.grid(True, linestyle='--', alpha=0.7)
            ax3.set_facecolor('white')
            ax3.legend(fontsize=10.5)
            ax3.tick_params(axis='x', labelsize=12)  # Set font size for x-axis tick labels
            ax3.tick_params(axis='y', labelsize=12)
        

