## Filename: fVParityPlot.py
## Date: 07/02/24
## Author: Itelcontar
## Description: Python3 module for creating Velocity Parity Plots, to determine quality of scaling of velocities by position, normalised as desired.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. noCellWidth - Integer number of cells to divide tank diameter linearly into. All voxels are cubic.
# 3. spuriousFilter - Integer, specifying minimum number of velocity data points in a pixel for the pixel to be included in the final plot. Removes 'spurious' velocities, with insufficient data for accurate averaging, or valid assumption of ergodicity.
# 4. norm - Determines how to 'normalise' velocities. Options are 'TankD' (Tank Diameter), 'TipSpeed' (Impeller Tip Velocity). Tip Speed is default.
        # Not case-sensitive. Both diameter, and tip speed calculated from data stored in globalData[].

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from FGlobalSettings import globalData

### FUNCTION ###

def fVParityPlot(exp,noCellWidth,spuriousFilter,norm):

    # Collecting data for both datasets
    posX1 = globalData[exp[0]]['posX']
    posY1 = globalData[exp[0]]['posY']
    posZ1 = globalData[exp[0]]['posZ']
    vX1 = globalData[exp[0]]['vX']
    vY1 = globalData[exp[0]]['vY']
    vZ1 = globalData[exp[0]]['vZ']
    TankD1 = globalData[exp[0]]['TankD']
    posX2 = globalData[exp[1]]['posX']
    posY2 = globalData[exp[1]]['posY']
    posZ2 = globalData[exp[1]]['posZ']
    vX2 = globalData[exp[1]]['vX']
    vY2 = globalData[exp[1]]['vY']
    vZ2 = globalData[exp[1]]['vZ']
    TankD2 = globalData[exp[1]]['TankD']
    name1 = globalData[exp[0]]['name']
    name2 = globalData[exp[1]]['name']
    
    # Next, find magnitudes, and normalise
    if norm.lower() == 'tipspeed':
        # Finding TipSpeeds
        radius1 = globalData[exp[0]]['TankD'] / 2
        RPM1 = globalData[exp[0]]['RPM']
        tipSpeed1 = (1/3) * radius1 * 2 * np.pi * RPM1 / 60 # (1/3) * 2piR * RPM / 60 => mm/s tip speed            
        radius2 = globalData[exp[1]]['TankD'] / 2
        RPM2 = globalData[exp[1]]['RPM']
        tipSpeed2 = (1/3) * radius2 * 2 * np.pi * RPM2 / 60 # (1/3) * 2piR * RPM / 60 => mm/s tip speed            

        # TipSpeed normalising
        vMag1 = np.divide(np.sqrt(vX1**2 + vY1**2 + vZ1**2), tipSpeed1)
        vMag2 = np.divide(np.sqrt(vX2**2 + vY2**2 + vZ2**2), tipSpeed2)
        units = "($v/V_{tip}$)"
    else:
        # TankD normalising
        vMag1 = np.divide(np.sqrt(vX1**2 + vY1**2 + vZ1**2), TankD1)
        vMag2 = np.divide(np.sqrt(vX2**2 + vY2**2 + vZ2**2), TankD2)
        units = "($TankD/s$)"
        
    # Finding a mutually beneficial 3D grid shape
    nx = noCellWidth
    cellSize = [((np.nanmax(posX1)+1)/noCellWidth), (np.nanmax(posX2)+1)/noCellWidth]
    ny1 = int(np.ceil(np.nanmax(posY1) / cellSize[0]))
    ny2 = int(np.ceil(np.nanmax(posY2) / cellSize[1]))
    nz1 = int(np.ceil(np.nanmax(posZ1) / cellSize[0]))
    nz2 = int(np.ceil(np.nanmax(posZ2) / cellSize[1]))
    ny = max(ny1,ny2)
    nz = max(nz1,nz2)

    # Defining function to create velocity scalar fields
    def calculateVelocityGrids(posX, posY, posZ, vMag, cellSize, nx, ny, nz):
        # Creating a rank 3 tensor for velocity magnitudes
        vGrid = np.zeros((ny, nx, nz))
        nGrid = np.zeros((ny, nx, nz))

        # Now run through the velocities, and store at the correct positions in the tensor
        for i in range(0,len(vMag)):
            x, y, z, v = posX[i], posY[i], posZ[i], vMag[i] # Aliasing for simplicity

            # Exclude data points with NaN values
            if not np.isnan(x) and not np.isnan(y) and not np.isnan(z) and not np.isnan(v):
                ix = int(x / cellSize)
                iy = int(y / cellSize)
                iz = int(z / cellSize)
                
                nGrid[iy, ix, iz] += 1  # Add one to the count
                vGrid[iy, ix, iz] += v  # Adding the velocity to the velocity count

        # To remove spurious velocities and carry out average
        vGrid[nGrid <= spuriousFilter] = np.nan
        vGrid = np.nan_to_num(vGrid / nGrid)  # Averaging
        
        return vGrid
        # END of General Function

    # Calculate velocity fields for both sets of data
    vGrid1 = calculateVelocityGrids(posX1, posY1, posZ1, vMag1, cellSize[0], nx, ny, nz)
    vGrid2 = calculateVelocityGrids(posX2, posY2, posZ2, vMag2, cellSize[1], nx, ny, nz)

    # Now we flatten the rank 3 tensors to column vectors
    vList1 = vGrid1.reshape(-1,1)
    vList2 = vGrid2.reshape(-1,1)

    # Dealign with NaNs
        # Check for NaN values in either vector
    mask = (np.isnan(vList1) | (vList1 < 0.01)) | (np.isnan(vList2) | (vList2 < 0.01))
    vList1 = np.array(vList1)[~mask]
    vList2 = np.array(vList2)[~mask]


    ### PLOTTING ###
    fig, ax = plt.subplots(figsize=(6, 6), tight_layout=True)

    ax.scatter(vList2,vList1, color='black', label = "Velocity Magnitudes", s = 5)
    ax.plot([min(vList2), max(vList2)], [min(vList2), max(vList2)], linestyle='--', color='red', linewidth=2, label = 'Parity')
    ax.set(ylabel=name1 + " Velocities " + units, xlabel=name2 + " Velocities " + units)
    ax.set_aspect(1)
    ax.legend(fontsize=10.5)

    # Statistical Measures
    rSquared = r2_score(vList1, vList2)
    fig.suptitle("Velocity Parity Plot - Cubic Grid Comparison " + str(noCellWidth) +"^3 Cells - R^2 = " + str(round(rSquared,2)),fontsize=14)

    return rSquared
