## Filename: FGlobalSettings.py
## Date: 22/12/23
## Author: Itelcontar
## Description: Constains all Global Definitions for Dictionary Data Structure. Essential for other PEPT scrips. Creates empty datastructure for data insertion.
        # MANUALLY vary datastructure length, as required for number of experiments.
        
### PARAMETERS ###
# 1. time - From CSV. List of particle detection times.
# 2. posX - From CSV. List of particle X positions.
# 3. posY - From CSV. List of particle Y positions.
# 4. posZ - From CSV. List of particle Z positions.
# 5. vX - From CSV. List of particle X velocities.
# 6. vY - From CSV. List of particle Y velocities.
# 7. vZ - From CSV. List of particle Z velocities.
# 8. TankD - Tank Diameter, in mm.
# 9. RPM - Tank RPM.
# 10. name - Name of run. This will appear on legends for plots, etc. 

### FUNCTION ###

def initialisation():
    global globalData
    globalData = {
        1: {'time': None, 'posX': None, 'posY': None, 'posZ': None, 'vX' : None, 'vY' : None, 'vZ' : None, 'TankD' : None, 'RPM' : None, 'name' : None},
        2: {'time': None, 'posX': None, 'posY': None, 'posZ': None, 'vX' : None, 'vY' : None, 'vZ' : None, 'TankD' : None, 'RPM' : None, 'name' : None},
        3: {'time': None, 'posX': None, 'posY': None, 'posZ': None, 'vX' : None, 'vY' : None, 'vZ' : None, 'TankD' : None, 'RPM' : None, 'name' : None},
        4: {'time': None, 'posX': None, 'posY': None, 'posZ': None, 'vX' : None, 'vY' : None, 'vZ' : None, 'TankD' : None, 'RPM' : None, 'name' : None},
        5: {'time': None, 'posX': None, 'posY': None, 'posZ': None, 'vX' : None, 'vY' : None, 'vZ' : None, 'TankD' : None, 'RPM' : None, 'name' : None},
        6: {'time': None, 'posX': None, 'posY': None, 'posZ': None, 'vX' : None, 'vY' : None, 'vZ' : None, 'TankD' : None, 'RPM' : None, 'name' : None},
        7: {'time': None, 'posX': None, 'posY': None, 'posZ': None, 'vX' : None, 'vY' : None, 'vZ' : None, 'TankD' : None, 'RPM' : None, 'name' : None},
        8: {'time': None, 'posX': None, 'posY': None, 'posZ': None, 'vX' : None, 'vY' : None, 'vZ' : None, 'TankD' : None, 'RPM' : None, 'name' : None},    
        9: {'time': None, 'posX': None, 'posY': None, 'posZ': None, 'vX' : None, 'vY' : None, 'vZ' : None, 'TankD' : None, 'RPM' : None, 'name' : None},
        10: {'time': None, 'posX': None, 'posY': None, 'posZ': None, 'vX' : None, 'vY' : None, 'vZ' : None, 'TankD' : None, 'RPM' : None, 'name' : None}
        }
