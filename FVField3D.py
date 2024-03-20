## Filename: FVField3D.py
## Date: 22/12/23
## Author: Itelcontar
## Description: Python3 module for creating 3D Velocity Vector Field Plots.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.
# 2. cellSize - In mm, single dimension size for cubic voxels.
# 3. spuriousFilter - Integer, specifying minimum number of velocity data points in a voxel for the voxel to be included in the final plot. Removes 'spurious' velocities, with insufficient data for accurate averaging, or valid assumption of ergodicity.
# 4. plots - Boolean numpy array for generating various plots. [MatPlotLib Quiver, mayAVI Quiver, MatPlotLib high-velocity voxel].     
# 5. animate - Similar to 'plots', except for animations. [mayAVI Vector field rotating view]
# NOTE - Recommend only creating one MAYAVI plot at once. 

### INTERNAL VARIABLES ###
# These should not need changing, but are a good starting point for troubleshooting.
AR = 'equal' # General purpose aspect ratio definition.
frames = 450 # Number of frames for an animation. 
framesFolder = 'MAYAVI Frames' # Name of folder for png animated frame storage (located in same location as this file).
vVmax = 0.4 # Fraction of maximum tank velocity necessary for a voxel to be included in high velocity voxel plots. Default 0.64.
onlyY = False # Manual variable to only consider y velocities in voxelplots.

import numpy as np
import matplotlib.pyplot as plt
from FGlobalSettings import globalData
from mpl_toolkits.mplot3d import Axes3D
from mayavi import mlab
import os
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel
from scipy.ndimage.filters import laplace # For divergence of velocity fields
from scipy.ndimage.filters import convolve # For vorticity

### FUNCTION ###

def fVField3D(exp, cellSize, spuriousFilter, plots, animate):

    # Running through 'exp'
    for i in range(0,len(exp)):

        # Finding TipSpeeds
        radius = globalData[exp[i]]['TankD'] / 2
        RPM = globalData[exp[i]]['RPM']
        tipSpeed = (1/3) * radius * 2 * np.pi * RPM / 60 # (1/3) * 2piR * RPM / 60 => mm/s tip speed            

        # Gathering data - all data is normalised by default
        TankD = globalData[exp[i]]['TankD']
        time = globalData[exp[i]]['time']
        posX = np.divide(globalData[exp[i]]['posX'],TankD)
        posY = np.divide(globalData[exp[i]]['posY'],TankD)
        posZ = np.divide(globalData[exp[i]]['posZ'],TankD)
        vX = np.divide(globalData[exp[i]]['vX'],tipSpeed)
        vY = np.divide(globalData[exp[i]]['vY'],tipSpeed)
        vZ = np.divide(globalData[exp[i]]['vZ'],tipSpeed)
        name = globalData[exp[i]]['name']

        # No cells in Grid.
        nx = int(np.ceil(np.nanmax(posX)*TankD / (cellSize)))
        ny = int(np.ceil(np.nanmax(posY)*TankD / (cellSize)))
        nz = int(np.ceil(np.nanmax(posZ)*TankD / (cellSize)))

        # Creating Grids - All same
        vxgrid = np.zeros((nx, ny, nz))
        vygrid = np.zeros((nx, ny, nz))
        vzgrid = np.zeros((nx, ny, nz))
        vMagGrid = np.zeros((nx, ny, nz))
        ngrid = np.zeros((nx, ny, nz))

        # For every velocity, identify location and record components 
        for i in range(len(vX)):
            x, y, z, vx, vy, vz = posX[i], posY[i], posZ[i], vX[i], vY[i], vZ[i]

            # Exclude data points with NaN values
            if not np.isnan(x) and not np.isnan(y) and not np.isnan(vx) and not np.isnan(vy):
                ix = int(x *TankD/ (cellSize))
                iy = int(y *TankD/ (cellSize))
                iz = int(z *TankD/ (cellSize))

                ngrid[ix, iy, iz] += 1  # Add one to the count
                vxgrid[ix, iy, iz] += vx  # Adding the velocity to the velocity count
                vygrid[ix, iy, iz] += vy
                vzgrid[ix, iy, iz] += vz
                vMagGrid[ix, iy, iz] += np.sqrt(vx**2 + vy**2 + vz**2)

        # To remove spurious velocities
        vxgrid[ngrid <= spuriousFilter] = np.nan
        vygrid[ngrid <= spuriousFilter] = np.nan
        vzgrid[ngrid <= spuriousFilter] = np.nan
        vMagGrid[ngrid <= spuriousFilter] = np.nan

        # Averaging velocities
        vxgrid = np.nan_to_num(vxgrid / ngrid)
        vygrid = np.nan_to_num(vygrid / ngrid)
        vzgrid = np.nan_to_num(vzgrid / ngrid)
        vMagGrid = np.nan_to_num(vMagGrid / ngrid)

        X, Y, Z = np.meshgrid(np.linspace(np.nanmin(posX), np.nanmax(posX), nx),
                              np.linspace(np.nanmin(posY), np.nanmax(posY), ny),
                              np.linspace(np.nanmin(posZ), np.nanmax(posZ), nz),
                              indexing='ij') # Specify 'ij', because PEPT 'Y' axis, is mathematical 'Z' axis
        # Saving Clean Grids
        vxgridClean = vxgrid
        vygridClean = vygrid
        vzgridClean = vzgrid
        vMagGridClean = vMagGrid
        
        ### PLOTTING ###

        if plots[0] == True:
            # Create 3D Vector Field Plot - MatPlotLib
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
                # N.B PEPT 'Y' axis, is mathematical 'Z' axis
            ax.quiver(X,Z,Y, vxgrid, vzgrid, vygrid, color='b', length=0.05, normalize=True, arrow_length_ratio=0.3, linewidth=0.3)

            # Plot Aesthetics
            ax.set_xlabel('X Position')
            ax.set_ylabel('Z Position')
            ax.set_zlabel('Y Position')
            ax.set_title("3D Vector Field Plot - " + name)
            ax.set_aspect(AR)

        if plots[1] == True:
            # Create Mayavi Plots
            fig1 = mlab.figure(size=(1500, 1200))
            
            # Removing all but the highest magnitude velocities
            vxgrid[vMagGrid <= vVmax*np.nanmax(vMagGrid)] = np.nan
            vygrid[vMagGrid <= vVmax*np.nanmax(vMagGrid)] = np.nan
            vzgrid[vMagGrid <= vVmax*np.nanmax(vMagGrid)] = np.nan

            # N.B PEPT 'Y' axis, is mathematical 'Z' axis
            f = mlab.quiver3d(X,Z,Y,vxgrid, vzgrid, vygrid, colormap = 'plasma')
            
            # Plot Aesthetics
            title = mlab.text(0.6, 0.9, name, width=0.3)
            #scene = mlab.gcf().scene
            #scene.background = (1,1,1)

            # Make a frames folder, if one doesn't already exist
            if not os.path.exists(framesFolder):
                os.makedirs(framesFolder)
            if not os.path.exists(framesFolder + "/Plot1/"):
                os.makedirs(framesFolder + "/Plot1/")
                
            if animate[0] == True:
                # Animating Vector Field Quiver, with rotation
                for n in range(frames):
                    # Define the camera position and orientation for each frame
                    azim = 360 * (n/frames)  # Adjust the angle as a function of time
                    mlab.view(azimuth=azim, elevation=-90 + 30*np.sin(azim*0.02), distance=2.5)

                    # Save each frame as an image
                    mlab.savefig(framesFolder + "/Plot1/" + f'frame_{n:03d}.png')

        # Removing all but the highest magnitude velocities
        vMagGrid[vMagGrid <= vVmax*np.nanmax(vMagGrid)] = np.nan

        if plots[2] == True:
            # VoxelPlot using Matplotlib:
            Voxel1MagGrid = np.logical_not(np.isnan(vMagGrid)) # Highest Velocities, ready to voxel plot

            if onlyY == True:
                 # Removing all but the highest magnitude velocities
                vygrid[vygrid <= vVmax*np.nanmax(vygrid)] = np.nan
                Voxel1MagGrid = np.logical_not(np.isnan(vygrid))
                #V = vygrid.flatten()
                name = name + " - Vy Only"

            # Plotting
            fig2 = plt.figure(figsize=(12, 12))
            ax2 = fig2.add_subplot(111, projection='3d')
            ax2.voxels(np.transpose(Voxel1MagGrid, (0, 2, 1))) # N.B PEPT 'Y' axis, is mathematical 'Z' axis

           # Plot Aesthetics
            ax2.set_title("3D Voxel Plot - Highest Velocities - " + name)
            ax2.set_aspect(AR)
            ax2.set_xlabel('X Position')
            ax2.set_ylabel('Z Position')
            ax2.set_zlabel('Y Position')
            ax2.set_aspect(AR)

        if plots[3] == True:
            # MayAVI Voxelplot
            x, z, y = np.indices(vMagGrid.shape)
            V = vMagGrid.flatten()

            if onlyY == True:
                 # Removing all but the highest magnitude velocities
                vygrid[vygrid <= vVmax*np.nanmax(vygrid)] = np.nan
                V = vygrid.flatten()
                name = name + " - Vy Only"
            
            # Plotting
            fig3 = mlab.figure(size=(1500, 1200))
            points = mlab.points3d(x.flatten(), y.flatten(), z.flatten(), V, mode='cube', colormap='viridis')
            mlab.colorbar(points, title='($v/V_{tip}$)')

            if animate[1] == True:
                # Animating Voxel Plot, with rotation
                
                # Make a frames folder, if one doesn't already exist
                if not os.path.exists(framesFolder):
                    os.makedirs(framesFolder)
                if not os.path.exists(framesFolder + "/Plot2/"):
                    os.makedirs(framesFolder + "/Plot2/")
                
                for n in range(frames):
                    # Define the camera position and orientation for each frame
                    azim = 360 * (n/frames)  # Adjust the angle as a function of time
                    mlab.view(azimuth=azim, elevation=50, distance=150)

                    #Move the camera center up a tad
                    mlab.move(0, 0, -5)
                    
                    # Save each frame as an image
                    mlab.savefig(framesFolder + "/Plot2/" + f'frame_{n:03d}.png')

        if plots[4] == True:
            # MayAVI Voxelplot Alternative
            VGrid = vMagGrid

            if onlyY == True:
                 # Removing all but the highest magnitude velocities
                vygrid[vygrid <= vVmax*np.nanmax(vygrid)] = np.nan
                VGrid = vygrid
                name = name + " - Vy Only"
                
                VGrid = rgrid

            # Generate data using nested loops
            x, y, z, V = [], [], [], []

            for i in range(VGrid.shape[0]):
                for j in range(VGrid.shape[1]):
                    for k in range(VGrid.shape[2]):
                        x.append(i)
                        y.append(j)
                        z.append(k)
                        V.append(VGrid[i, j, k])

            # Convert lists to numpy arrays
            x, y, z, V = np.array(x), np.array(y), np.array(z), np.array(V)
            
            # Mask NaN values
            validIndices = ~np.isnan(V)
            x = x[validIndices]
            y = y[validIndices]
            z = z[validIndices]
            V = V[validIndices]

            fig4 = mlab.figure(size=(1500, 1200))
            points = mlab.points3d(x,z,y,V, mode='cube', scale_mode = 'none', scale_factor = 1.5)

            points.glyph.color_mode = 'color_by_scalar'
            cbar = mlab.colorbar(points, title='($v/V_{tip}$)')
            cbar.scalar_bar.width = 0.4  # Adjust the width as needed
            cbar.scalar_bar.height = 0.05  # Adjust the height as needed
            cbar.label_text_property.font_size = 5
            cbar.label_text_property.color = (0, 0, 0) # Text to black
            cbar.title_text_property.color = (0, 0, 0)
            title = mlab.text(0.6, 0.9, name, width=0.3)

            scene = mlab.gcf().scene
            scene.background = (1,1,1)
            
            if animate[2] == True:
                # Animating Voxel Plot, with rotation
                
                # Make a frames folder, if one doesn't already exist
                if not os.path.exists(framesFolder):
                    os.makedirs(framesFolder)
                if not os.path.exists(framesFolder + "/Plot2/"):
                    os.makedirs(framesFolder + "/Plot2/")
                
                for n in range(frames):
                    # Define the camera position and orientation for each frame
                    azim = 360 * (n/frames)  # Adjust the angle as a function of time
                    mlab.view(azimuth=azim, elevation=-80 + 13*np.sin(azim*0.05), distance=150)

                    #Move the camera center up a tad
                    mlab.move(0, 0, -5)
                    
                    # Save each frame as an image
                    mlab.savefig(framesFolder + "/Plot2/" + f'frame_{n:03d}.png')

        # High Vorticity Plot
        if plots[5] == True:
            # Collecting Data
            vxGrid5 = vxgridClean
            vyGrid5 = vygridClean
            vzGrid5 = vzgridClean

            # Correct computation of vorticity components
            dvz_dy, dvz_dx, dvz_dz = np.gradient(vzGrid5)
            dwy_dx, dwy_dy, dwy_dz = np.gradient(vyGrid5)
            dvx_dz, dvx_dy, dvx_dx = np.gradient(vxGrid5)

            vorticity_x = dwy_dz - dvz_dy
            vorticity_y = dvx_dx - dwy_dx
            vorticity_z = dvz_dy - dvx_dy

            # Handle zeros in the denominator
            epsilon = np.finfo(float).eps
            dvx_dx[dvx_dx == 0] = epsilon
            dwy_dy[dwy_dy == 0] = epsilon
            dvz_dz[dvz_dz == 0] = epsilon
            vorticityMagnitude = np.sqrt(vorticity_x**2 + vorticity_z**2)
            VGrid = vorticityMagnitude

            # Thresholding
            VGrid[VGrid < 0.4*np.nanmax(VGrid)] = np.nan
            VGrid[VGrid > 0.9*np.nanmax(VGrid)] = np.nan

            # Create Mayavi Plots
            fig5 = mlab.figure(size=(1500, 1200))

            # N.B PEPT 'Y' axis, is mathematical 'Z' axis
            # Generate data using nested loops
            x, y, z, V = [], [], [], []

            for i in range(VGrid.shape[0]):
                for j in range(VGrid.shape[1]):
                    for k in range(VGrid.shape[2]):
                        x.append(i)
                        y.append(j)
                        z.append(k)
                        V.append(VGrid[i, j, k])

            # Convert lists to numpy arrays
            x, y, z, V = np.array(x), np.array(y), np.array(z), np.array(V)
            
            # Mask NaN values
            validIndices = ~np.isnan(V)
            x = x[validIndices]
            y = y[validIndices]
            z = z[validIndices]
            V = V[validIndices]

            fig4 = mlab.figure(size=(1500, 1200))
            points = mlab.points3d(x,z,y,V, mode='cube', scale_mode = 'none', scale_factor = 1)

            points.glyph.color_mode = 'color_by_scalar'
            cbar = mlab.colorbar(points, title='($v/V_{tip}$)')
            cbar.scalar_bar.width = 0.4  # Adjust the width as needed
            cbar.scalar_bar.height = 0.05  # Adjust the height as needed
            cbar.label_text_property.font_size = 5
            cbar.label_text_property.color = (0, 0, 0) # Text to black
            cbar.title_text_property.color = (0, 0, 0)
            title = mlab.text(0.6, 0.9, name, width=0.3)

            scene = mlab.gcf().scene
            scene.background = (1,1,1)
            # Plot Aesthetics
            title = mlab.text(0.6, 0.9, name, width=0.3)

            # Make a frames folder, if one doesn't already exist
            if not os.path.exists(framesFolder):
                os.makedirs(framesFolder)
            if not os.path.exists(framesFolder + "/Plot1/"):
                os.makedirs(framesFolder + "/Plot1/")
                
            if animate[3] == True:
                # Animating Vector Field Quiver, with rotation
                for n in range(frames):
                    # Define the camera position and orientation for each frame
                    azim = 360 * (n/frames)  # Adjust the angle as a function of time
                    mlab.view(azimuth=azim, elevation=-90 + 10*np.sin(azim*0.05), distance=2.5)

                    # Save each frame as an image
                    mlab.savefig(framesFolder + "/Plot1/" + f'frame_{n:03d}.png')
            

            
