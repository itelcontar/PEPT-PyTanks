## Filename: FCleanPlot.py
## Date: 30/12/23
## Author: Itelcontar
## Description: Python3 module for creating plots, highlighting NaN data, removed from PEPT data during filtering/cleaning process.

### PARAMETERS ###
# 1. exp - Numpy array of integers corresponding to experiments stored in globalData[] dictionary. See 'FDataPull' for further information.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from FGlobalSettings import globalData

### FUNCTION ###

def fCleanPlot(exp):
    for i in range(0,len(exp)):

        # Gathering data
        time = globalData[exp[i]]['time']
        posX = globalData[exp[i]]['posX']
        posY = globalData[exp[i]]['posY']
        posZ = globalData[exp[i]]['posZ']
        vX = globalData[exp[i]]['vX']
        vY = globalData[exp[i]]['vY']
        vZ = globalData[exp[i]]['vZ']
        name = globalData[exp[i]]['name']

        # Masking NaNs
        TotalMask = np.full(len(posX),True)
        TotalMask = np.logical_not(np.isnan(vX)) #Contribs from velocity calc boundaries (gets messy with Data Cleaning)
        FinalMask = np.logical_not(np.isnan(posX))
        TotalMask = np.logical_and(np.logical_not(np.isnan(posX)),FinalMask) # Contribs from Position Data

        # Identify the ranges of excluded data
        excluded_ranges = []
        start_idx = None
        for i, is_valid in enumerate(TotalMask):
            if is_valid and start_idx is not None:
                excluded_ranges.append((start_idx, i - 1))
                start_idx = None
            elif not is_valid and start_idx is None:
                start_idx = i

        # If the last range extends to the end of the data
        if start_idx is not None:
            excluded_ranges.append((start_idx, len(time) - 1))
        print("\nIndex ranges excluded from data: \n" + str(excluded_ranges)+ "\n")

        # Plotting using MatPlotLib
        figs, axs = plt.subplots(3,1,figsize = (20,4), sharey=True, sharex= True, tight_layout = True)
        line_posX, = axs[0].plot(time/3600,posX, label='1D Particle Position Data')
        axs[0].set(ylabel = 'x $(mm)$')
        axs[1].plot(time/3600,posY)
        axs[1].set(ylabel = 'y $(mm)$')
        axs[2].plot(time/3600,posZ)
        axs[2].set(xlabel='Run Time $(hr)$', ylabel = 'z $(mm)$')
        figs.suptitle("Particle Position with Cleaned Data: " + name, fontsize = 14)

        # Adding red Rectangles to show removed data
        for start, end in excluded_ranges:
            for ax in axs:
                ax.add_patch(patches.Rectangle(
                    (time[start] / 3600, ax.get_ylim()[0]),
                    (time[end] - time[start]) / 3600,
                    ax.get_ylim()[1] - ax.get_ylim()[0],
                    linewidth=0, edgecolor='none', facecolor='red', alpha=0.5
                    ))

        # Plot aesthetics
        for ax in axs:
            ax.grid(True, linestyle='--', alpha=0.7)
            ax.set_facecolor('#f0f0f0')

        legend_rectangles = [patches.Rectangle((0, 0), 1, 1, facecolor='red', alpha=0.5, edgecolor='none', label='Excluded Data')]
        axs[0].legend(handles=[line_posX] + legend_rectangles, loc='upper right') 

        
