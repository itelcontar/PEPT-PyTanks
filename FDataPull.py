## Filename: FDataPull.py
## Date: 22/12/23
## Author: Itelcontar
## Description: Imports data to dictionary structure into a HMI master file.

### PARAMETERS ###
# 1. path - Path to csv file of PEPT data, not including filename.
# 2. file - Filename of csv file, with lagrangian data in columns time, posX, posY, posZ. E.g. 'run1.csv' 
# 3. fileNumber - Number to enter into dictionary structure. Referenced as 'epx' in other scripts retrieving data. 
# 4. TankD - Tank Diameter, in mm.
# 5. RPM - Tank RPM.
# 6. name - Name of run. This will appear on legends for plots, etc. 

### FUNCTION ###

import csv
import numpy as np
from FGlobalSettings import globalData

def fDataPull(path,file,fileNumber,TankD,RPM,name):

    # Collecting CSV Data
    data = np.genfromtxt(path + file, delimiter=',')

    # In the event of no pre-allocated space in FGlobalSettings
    if fileNumber not in globalData:
        globalData[fileNumber] = {}

    # Assigning data to names in dictionary data-structure
    globalData[fileNumber]['time'] = data[:, 0]
    globalData[fileNumber]['posX'] = data[:, 1]
    globalData[fileNumber]['posY'] = data[:, 2]
    globalData[fileNumber]['posZ'] = data[:, 3]
    globalData[fileNumber]['vX'] = data[:, 4]
    globalData[fileNumber]['vY'] = data[:, 5]
    globalData[fileNumber]['vZ'] = data[:, 6]
    globalData[fileNumber]['TankD'] = TankD
    globalData[fileNumber]['RPM'] = RPM
    globalData[fileNumber]['name'] = name
    
