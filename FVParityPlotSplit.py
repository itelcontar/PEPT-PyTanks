## Filename: FVParityPlotSplit.py
## Date: 07/02/24
## Author: Itelcontar
## Description: Python3 module for creating impeller plane split velocity Parity Plots, to determine quality of scaling of KE by position, normalised as desired.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Only two experiments please. Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. noCellWidth - Integer number of cells to divide tank diameter linearly into. All voxels are cubic.
# 3. spuriousFilter - Integer, specifying minimum number of velocity data points in a pixel for the pixel to be included in the final plot. Removes 'spurious' velocities, with insufficient data for accurate averaging, or valid assumption of ergodicity.
# 4. DA - Depth Average? Deternime whether pizels or voxels are used. Voxels may be more true to reality, but pixels lack the noise present in voxel data, due to avraging. Choose carefully.
# 5. int00 - Boolean. Indicates whether best fit line should pass through origin.
# 6. norm - Determines how to 'normalise' velocities. Options are 'TankD' (Tank Diameter), 'TipSpeed' (Impeller Tip Velocity). Tip Speed is default.
        # Not case-sensitive. Both diameter, and tip speed calculated from data stored in globalData[].
# 7. plots - Boolena indicating which plots to create. [Parity Plot, Scalar Fields for Checking Comparison]

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from FGlobalSettings import globalData

### FUNCTION ###

def fVParityPlotSplit(exp,noCellWidth,spuriousFilter,DA,int00,norm,plots):

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
        tipKE1 = tipSpeed1 * tipSpeed1
        tipKE2 = tipSpeed2 * tipSpeed2
        
        # TipSpeed normalising
        vMag1 = np.divide(np.sqrt(vX1**2 + vY1**2 + vZ1**2), tipSpeed1)
        vMag2 = np.divide(np.sqrt(vX2**2 + vY2**2 + vZ2**2), tipSpeed2)
        units = "($v/{V_{tip}}$)"
        
        #units = "($v^{2}/{V_{tip}}^{2}$)"
        print("Running Tip Speed Normalisation")
    else:
        # TankD normalising
        vMag1 = np.divide((vX1**2 + vY1**2 + vZ1**2), TankD1)
        vMag2 = np.divide((vX2**2 + vY2**2 + vZ2**2), TankD2)
        units = "($TankD/s$)"
        
    # Ensuring no out-of-volume data exists that would disrup the plot:
    posX1 = np.where(abs(posX1-np.nanmean(posX1)) > TankD1/2, np.nan, posX1)
    posZ1 = np.where(abs(posZ1-np.nanmean(posZ1)) > TankD1/2, np.nan, posZ1)
    posX2 = np.where(abs(posX2-np.nanmean(posX2)) > TankD2/2, np.nan, posX2)
    posZ2 = np.where(abs(posZ2-np.nanmean(posZ2)) > TankD2/2, np.nan, posZ2)
       
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
    
    # Get the shape of the tensor
    shape = vGrid1.shape

    # Define the split index along the x-plane
    splitIndexY = int(shape[0] / 3)  # Change the denominator as needed

    # Split the tensor along the Y impeller-plane
    vGrid1Top = vGrid1[:splitIndexY, :, :]
    vGrid2Top = vGrid2[:splitIndexY, :, :]

    vGrid1Bottom = vGrid1[splitIndexY:, :, :]
    vGrid2Bottom = vGrid2[splitIndexY:, :, :]

    # Depth Averaging?
    if DA == True:
        vGrid1Top = np.nanmean(vGrid1Top,axis = 2)
        vGrid2Top = np.nanmean(vGrid2Top,axis = 2)
        vGrid1Bottom = np.nanmean(vGrid1Bottom,axis = 2)
        vGrid2Bottom = np.nanmean(vGrid2Bottom,axis = 2)
        vGrid1 = np.nanmean(vGrid1,axis = 2)
        vGrid2 = np.nanmean(vGrid2,axis = 2)

    # Now flatten the split tensors
    vList1Top = vGrid1Top.reshape(-1, 1)
    vList2Top = vGrid2Top.reshape(-1, 1)

    vList1Bottom = vGrid1Bottom.reshape(-1, 1)
    vList2Bottom = vGrid2Bottom.reshape(-1, 1)

    # Continue with NaN removal
    maskTop = (np.isnan(vList1Top) | (vList1Top < 0.0001)) | (np.isnan(vList2Top) | (vList2Top < 0.0001))
    maskBottom = (np.isnan(vList1Bottom) | (vList1Bottom < 0.0001)) | (np.isnan(vList2Bottom) | (vList2Bottom < 0.0001))

    vList1Top = np.array(vList1Top)[~maskTop]
    vList2Top = np.array(vList2Top)[~maskTop]
    vList1Bottom = np.array(vList1Bottom)[~maskBottom]
    vList2Bottom = np.array(vList2Bottom)[~maskBottom]

    # Plot 1 - Parity Plot
    if plots[0] == True:
        # Plotting both set of data differently on parity plot.
        fig, ax = plt.subplots(figsize=(6, 6), tight_layout=True)

        ax.scatter(vList2Top, vList1Top, color='black', label="Above Impeller Plane", s=5)
        ax.scatter(vList2Bottom, vList1Bottom, color='blue', label="Below Impeller Plane", s=5)

        # Parity line & Aesthetics
        ax.plot([0, max(np.nanmax(vList2Top),np.nanmax(vList2Bottom))], [0, max(np.nanmax(vList2Top),np.nanmax(vList2Bottom))], linestyle='--', color='red',
                linewidth=2, label='Parity')

        ax.set_ylabel(name1 + " V " + units, fontsize=14)
        ax.set_xlabel(name2 + " V " + units, fontsize=14)
        ax.legend()
        ax.tick_params(axis='x', labelsize=12)  # Set font size for x-axis tick labels
        ax.tick_params(axis='y', labelsize=12)

        # PolyFit for Lines of Best Fit
            # Default Model
        if int00 != True:
            modelTop = np.polyfit(vList2Top,vList1Top,1)
            modelBottom = np.polyfit(vList2Bottom,vList1Bottom,1)
            vHatTop = np.polyval(modelTop,vList2Top)
            vHatBottom = np.polyval(modelBottom,vList2Bottom)
            print("Above Impeller Plane Model " + "V2 = " + str(modelTop[0]) + " V1  +  " + str(modelTop[1]))
            print("Below Impeller Plane Model " + "V2 = " + str(modelBottom[0]) + " V1  +  " + str(modelBottom[1]))
        else: # Model passing through 0,0
            modelTop = np.sum(vList1Top * vList2Top) / np.sum(vList2Top**2)
            modelBottom = np.sum(vList1Bottom * vList2Bottom) / np.sum(vList2Bottom**2)
            vHatTop = modelTop * vList2Top
            vHatBottom = modelBottom * vList2Bottom
            print("Above Impeller Plane Model " + "V2 = " + str(modelTop))
            print("Below Impeller Plane Model " + "V2 = " + str(modelBottom))
            
        # Plotting
        ax.plot(vList2Top, vHatTop, color='black', linestyle=':', linewidth=2, label='Best-fit Line Above Impeller Plane')
        ax.plot(vList2Bottom, vHatBottom, color='blue', linestyle=':', linewidth=2, label='Best-fit Line Below Impeller Plane')

        # Statistical Measures for both portions
        rSquaredTop = r2_score(vList1Top, vList2Top)
        rSquaredBottom = r2_score(vList1Bottom, vList2Bottom)
        
        fig.suptitle("V Parity Plot " + str(noCellWidth) + "R^2 Above = " +
                     str(round(rSquaredTop, 2)) + ", R^2 Below = " + str(round(rSquaredBottom, 2)), fontsize=14)
        # R^2 Printing:
            #R62 for total parity
        set1 = np.append(vList1Top,vList1Bottom)
        set2 = np.append(vList2Top,vList2Bottom)
        r2All = r2_score(set1,set2)
        print("R^2t to Parity: " + str(r2All))

        
    # Scalar field velocity plot for Checking correlation:
    if plots[1] == True:
        # If not already 2D, make it so
        if DA == False:
            imGrid1 = np.nanmean(vGrid1,axis = 2)
            imGrid2 = np.nanmean(vGrid2,axis = 2)
        else:
            imGrid1 = vGrid1
            imGrid2 = vGrid2
            
        # Define figures
        fig1, ax1 = plt.subplots(figsize=(6, 6), tight_layout=True)
        ax1.imshow(imGrid1, origin = 'lower', extent=[np.nanmin(posX1),np.nanmax(posX1),np.nanmin(posY1),np.nanmax(posY1)], cmap='plasma')
        fig2, ax2 = plt.subplots(figsize=(6, 6), tight_layout=True)
        ax2.imshow(imGrid2, origin = 'lower', extent=[np.nanmin(posX2),np.nanmax(posX2),np.nanmin(posY2),np.nanmax(posY2)], cmap='plasma')
       
    # Optionally return both r-squared values
    return rSquaredTop, rSquaredBottom
