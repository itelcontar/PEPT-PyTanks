## Filename: PEPTPostprocessing.py
## Date: 16/01/24
## Author: Itelcontar
## Description: Imports cleaned data, and pulls plots from function files. Control panel for analysis.
## N.B. Thanks to Windows-Yule et. al for guidance in 'A Comprehensive Guide to PEPT'.

### INTRODUCTION ###
# This py file is an example of a master file, that demonstrates how to interact with the many different python
# scripts developed, for visualising PEPT data. Each of the scripts available is annotated, and preceeded with
# variable definitions, which are not repeated here.

# Input files must be csv, space-separated (see below to change to comma), containing columns of filtered PEPT data
# from left-to-right: 'time', 'X position', 'Y position', 'Z position', 'X velocity', 'Y velocity', 'Z velocity'.
# Data must be ALREADY filtered to remove positions outside the vessel volume. Another script has been provided, which
# includies this function. Data filtered out in this way should be replaced with NaNs. Time data should not be removed
# at all. 

### Setting Up ###
##
##

import numpy as np                  # Base requirement
import time as t                    # Useful, not essential
import matplotlib.pyplot as plt     # Base requirement
from mayavi import mlab             # Required for 3D plotting only

print("Filename: PEPTPostprocessing.py \nAuthor: Itelcontar")
print("\nScript Run starting ... ... ...\n")
start_time = t.time()

##
### I: Data Collection ###
##
##
##

print("Starting data collection... ")

# CSV File Names - for data extraction 
path = "InsertPathHere"
file1 = "1_S_Base_Ungassed.csv"                 # Small Ungassed
file2 = "2_S_Base_3VVM.csv"                     # Small Gassed
file3 = "3_M_Ungassed.csv"                      # Medium Ungassed
file4 = "4_M_VSG.csv"                           # Medium VSG VPD Scaling
file5 = "5_M_VVM.csv"                           # Medium VVM VPD Scaling
file8 = "8_L_VSG.csv"                           # Large VSG VPD Scaling
file9 = "9_L_Ungassed.csv"

names = ["T = 200mm, Ungassed", "T = 200mm, Gassed", "T = 288mm, Ungassed P/V", "T = 288 mm, VSG & P/V", "T = 288 mm, VVM & P/V", "DEV. Large UG 1","DEV. Large UG 2", "T = 450 mm, VSG & P/V", "T= 450 mm, Ungassed P/V"]

# Dictionary of all File Results
from FGlobalSettings import initialisation
initialisation()
from FGlobalSettings import globalData

# Pulling in Data
from FDataPull import fDataPull # fdataPull(path_to_file, identification name, identification No., tank size (mm), impeller RPM, legend name)

fDataPull(path,file1,1,200,393,names[0])
fDataPull(path,file2,2,200,530,names[1])
fDataPull(path,file3,3,288,308,names[2])
fDataPull(path,file4,4,288,407,names[3])
fDataPull(path,file5,5,288,409,names[4])
fDataPull(path,file8,8,450,295,names[7])
fDataPull(path,file9,9,450,229,names[8])

print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()


##
### 1: Single Velocity Distirbution Plots
##
##
##

print("Starting distribution plotting...")

exp = []
bins = 100
plots = [False,True]

if exp != []:
    from FVDistribSingle import fVDistribSingle
    fVDistribSingle(exp, bins, '', 'vY', plots)
    fVDistribSingle(exp, bins, 'tipspeed', 'vX', plots)
    fVDistribSingle(exp, bins, 'tipspeed', 'v', plots)
    fVDistribSingle(exp, bins, 'tipspeed', 'vz', plots)

print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()
    
##
### 2: Residency Plots
##
##
##

print("Starting residency plotting...")

exp = [1]
cellSizes = [5,5,5]
angles = [1,2,3]


if exp != []:
    from FResid import fResid
    fResid(exp,cellSizes,angles)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()

##
### 3: Velocity Fields
##
##
##

print("Starting velocity field plotting...")

exp = []
cellSize = 9 * 288 / 200
angles = [1,2,3,1,2,3]
spuriousFilter = 10

if exp != []:
    from FVField import fVField
    fVField(exp, cellSize, angles, spuriousFilter)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()

##
### 4: Velocity Scalar Fields
##
##
##

print("Starting scalar field plotting...")

exp = []
cellSize = 3 
angles = [1,2,3,1,2,3]
spuriousFilter = 2

if exp != []:
    from FVFieldScal import fVFieldScal
    fVFieldScal(exp, cellSize, angles, spuriousFilter)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()

##
### 5: Triple Scalar Fields
##
##
##

print("Starting triple scalar field plotting...")

exp = []
cellSize = 2*450/200
spuriousFilter = 4
plots = [False,True,False] # [scalarV, scalarT, vectorV]

if exp != []:
    from FTripleScalarField import fTripleScalarField
    fTripleScalarField(exp, cellSize, spuriousFilter, plots)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()


##
### 6: Velocity Distributions
##
##
##

print("Starting velocity distribution plotting...")

exp = []
bins = 100
norm = '' # 'TankD' or 'TipSpeed' or ''
plots = [True,True] # 

from FVDistrib import fVDistrib
if exp != []:
    fVDistrib(exp,bins,norm,plots)
    fVDistrib(exp,bins,'TankD',plots)
    fVDistrib(exp,bins,'TipSpeed',plots)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()


##
### 7: Differential Velocity Fields
##
##
##

print("Starting velocity differential plotting...")

exp = [] # Only Two Inputs Please!!
noCellWidth = 31
plots = [True, False, True, False] # True/False [Superimposed VectorField, Error Field, Contour Field]
angle = 1
spuriousFilter = 20

from FVFieldDiff import fVFieldDiff
if exp != []:
    fVFieldDiff(exp, noCellWidth, 1, spuriousFilter, plots)
    fVFieldDiff(exp, noCellWidth, 2, spuriousFilter, plots)
    fVFieldDiff(exp, noCellWidth, 3, spuriousFilter, plots)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()

##
### 8: Filter Assessment Plot
##
##
##

print("Starting plot filter data plotting...")

exp = []

from FCleanPlot import fCleanPlot
if exp != []:
    fCleanPlot(exp)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()

##
### 9: 1D Residency Plot
##
##
##

print("Starting 1D residency plotting...")

exp = []
angle = "posY"
nPoints = 40
norm = True
if exp != []:
    from FResid1D import fResid1D
    fResid1D(exp,angle,norm,nPoints)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()

##
### 10: Distributed Distributions of Velocities
##
##
##

print("Starting velocity distribution distributions plotting...")

exp = []
angle = 1
segments = 10
velocitySpec = 'vX' # vX vY vZ or v
bins = 100
norm = False

if exp != []:
    from FVDistribSlice import fVDistribSlice
    fVDistribSlice(exp, angle, segments, velocitySpec, bins, norm)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()


##
### 11: Traj Plotting
##
##
##

print("Starting 2D trajectory plotting...")

exp = []
length = 30000 # no. Points
angle = [1,1,1,1,1,1]

if exp!= []:
    from FTraj import fTraj
    fTraj(exp, length, angle)
    fTraj(exp, length, angle)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()

##
### 12: Recirculation Time
##
##
##

print("Starting recirculation plotting...")

exp = []
cellSize = [0.01]
plane = 0.3333
plots = [False,True,True]
mlt = 25
noLoops = 2
PlaneTolerance = 0.17

if exp != []:
    from FRCTime import fRCTime
    fRCTime(exp, cellSize, plane, plots, mlt, noLoops,PlaneTolerance)

##
### 13: Bubble Velocities
##
##
##

print("Starting bubble velocities plotting...")

exp = []
mjt = 0.3 # Minimum Jump Time (s)
mjd = 0.5 # Minimum Jump Distance (tankD)
plots = [True, True, True]

if exp!= []:
    from FJumpV import fJumpV
    fJumpV(exp,mjt,mjd,plots)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()

##
### 14: Up-Down Cut Velocities
##
##
##

print("Starting Up-Down Cut Velociites plotting...")

exp = []
cellSize = [1.4*4.5/2,1.4*4.5/2,1.4*2.9/2] # 2/3 for S/M
cutPoint = 0.333 # Plane at which we cut the plots
component = 'Vy' # 'vX', 'vZ', 'vY', 'V' or 'vXZ' for sqrt(vX**2 + vZ**2) - best patterns with 'vY' or 
spuriousFilter = 2

if exp!= []:
    from FVUpDown import fVUpDown
    fVUpDown(exp,cellSize,component,spuriousFilter,cutPoint)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()

##
### 15: Velocity Parity Plot
##
##
##

print("Starting Velocity Parity plotting...")

exp = [] # Only 2 inputs, for comparison.
noCellWidth = 15
spuriousFilter = 2
norm = 'TipSpeed' # Default 'TipSpeed', alternative 'TankD' 

if exp!= []:
    from FVParityPlot import fVParityPlot
    fVParityPlot(exp,noCellWidth,spuriousFilter,norm)
    fVParityPlot(exp,noCellWidth,spuriousFilter,'TankD')

# Location optimising
exp = [] # Only 2 inputs, for comparison.
noCellWidth = 14
spuriousFilter = 4
norm = 'Tipspeed' # Default 'TipSpeed', alternative 'TankD' 
offset = 0.2 # Fraction of tank diameter for largest offset extent. Should be <0.05.
nCR = 30 # Cube root of total optimisaion scenarios run. Run time grows fast!

if exp!= []:
    from FVParityPlot_LocOptimise import fVParityPlot_LocOptimise
    fVParityPlot_LocOptimise(exp,noCellWidth,spuriousFilter,offset,norm,nCR)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()


##
### 16: 3D Fields
##
##
##

print("Starting 3D Vector plotting...")

exp = []
cellSize = 6*2.88/2 
spuriousFilter = 2
plots = [False, True, False, False, False, False] # [MPL Quiver 3D, MAYAVI Quiver3D, MPL Voxels, MAYAVI Voxels, MAYAVI Voxels II]
animate = [True, False, False, False]             # [] 

if exp!= []:
    from FVField3D import fVField3D
    fVField3D(exp,cellSize,spuriousFilter,plots,animate)

print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()


##
### 17 : Cylin. Velocities, A Priori
##
##
##

print("Starting cylindrical velocity plotting...")

exp = []
jValue = 10
bins = 250
norm = 'TipSpeed' # Set to [], 'TankD' or 'TipSpeed'.

if exp!= []:
    from FVCylin import fVCylin
    fVCylin(exp,jValue,bins,norm,12)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()

##
### 18: JValue Investigate
##
##
##

print("Starting jValue Investigation...")

exp = []
jValueRange = [1]

if exp!= []:
    from FJValueInvestigate import fJValueInvestigate
    fJValueInvestigate(exp,jValueRange,[True, True])

print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()



##
### 19: Cartesian Velocity Regeneration
##
##
##

print("Starting A-Priori Velocity Calculation...")

exp = []
jValue = 1

if exp!= []:
    from FVAPriori import fVAPriori
    fVAPriori(exp,jValue,[True, True])

print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()


##
### 20 : Data Rate
##
##
##

exp = []

if exp!= []:
    from FDataRate import fDataRate
    fDataRate(exp)

##
### 21 : Split Parity Plots
##
##
##

print("Starting impeller plane split parity plots...")

exp = [] # Only 2 inputs, for comparison.
noCellWidth = 25
spuriousFilter = 10
norm = 'Tipspeed' # Default 'TipSpeed', alternative 'TankD'
DA = True # Depth Average? 
int00 = True
plots = [True,False]

if exp!= []:
    from FVParityPlotSplit import fVParityPlotSplit
    fVParityPlotSplit(exp,noCellWidth,spuriousFilter,DA,int00,norm,plots)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()


##
### 22: Radial Velocity Plots
##
##
##

print("Starting radial velocity plots...")

exp = []
component = 'V'

if exp!= []:
    from FVRad import fVRad
    fVRad(exp,component)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()


##
### 23: Acceleration PLots
##
##
##

print("Starting acceleration plots...")

exp = []
cellSize = 7
angles = [3,3]
spuriousFilter = 5
plots = [True,True,True]

if exp!= []:
    from FAField import fAField
    fAField(exp, cellSize, angles, spuriousFilter,plots)
print("...completed in " + str(round(t.time() - start_time)) + " s.")
start_time = t.time()


##
### 24: Show All Plots
##
##
##

plt.show()
mlab.show()
print("\n... ... ... script run is complete!")


