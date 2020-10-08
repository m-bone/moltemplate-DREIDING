# Import packages
import os
import sys
import math
import numpy as np
import pandas as pd
import xyzFunctions as xf

directory = sys.argv[1]
name = sys.argv[2]
box = sys.argv[3]
ff = sys.argv[4]

# Set working directory
os.chdir(directory)

# Post-process user input
filename = name + ".xyz"
molecule = name.upper()
outputFile = name + ".lt"
forcefield = ff + ".lt"
FORCEFIELD = ff.upper() # forcefield file must use uppercase function

# Set up dictionaries
atomDictionary = xf.atomDictionary

# Import xyz file using np and convert to pd dataframe
atomData = xf.xyzToDF(filename)

# Create a list of tuples of (int, list) that contain the bonder (bondData[x][0]) and bondees (bondData[x][1][y])
bondData = xf.bondCalculation(atomData)
# Save current element list
elementList = list(atomData['Element'])

# Set for LAMMPS 'full' style
dataStyle = " $mol:... @atom:"

# Add necessary columns to dataframe
atomData.insert(1, "Jargon", value="NaN")
atomData.insert(2, "Charge", value=0.0)
atomData["Jargon"] = dataStyle

# Update atom names with counter values
# Return list of atom names
atomName = xf.moltemplater(atomData, atomDictionary)

# Create readable list of bonds between atoms
with open(molecule + '.bondlist.txt', 'w') as readableBonds:

    for index, value1 in enumerate(bondData):
        # Atom X ...
        readableBonds.write("Atom: " + atomName[value1[0]] + "\n")
        for value2 in bondData[index][1]:
            # ... bonded with Atoms Y,Z, etc
            readableBonds.write("Bonded with: " + atomName[value2] + "\n")
        readableBonds.write("\n")

# Write main file
################################################################################
# Open file
with open(outputFile, 'w') as file:

    # Begin molecule
    file.write('import "' + forcefield + '" \n\n')
    file.write(molecule + " inherits " + FORCEFIELD + " {\n")

    # Write atom data
    file.write('\twrite("Data Atoms") {\n')
    np.savetxt(file, atomData, fmt='\t\t %s %s %s %f %f %f')
    file.write('\n\t}\n')

    # Write bond data
    file.write('\twrite("Data Bond List") {\n')

    # Preload counter and empty np array for storing already made bonds
    bondCounter = 0
    prevBond = np.empty((0,2), int)

    # Write bonds section using earlier defined bond data
    for index, value1 in enumerate(bondData):
        for value2 in bondData[index][1]:
            # Bond pairs that have already been will return a length of 1
            repeat = np.where((prevBond == (value1[0], value2)).all(axis=1))[0]
            if len(repeat) != 0:
                continue # Stops bonds being repeated
            file.write('\t\t$bond:' + str(bondCounter) + '\t$atom:' + str(atomName[value1[0]]) + '\t$atom:' + str(atomName[value2]) + '\n')
            bondCounter += 1
            prevBond = np.append(prevBond, np.array([[value2, value1[0]]]), axis = 0) # Flip order around to catch when p and q are reversed
    file.write('\t}\n')

    # Close molecule bracket and save file
    file.write('}')


# Write system.lt file
################################################################################
with open("system.lt", "w") as system:

    # General setup of atom and cubic cell
    system.write('import "' + outputFile + '"\n\n')
    system.write("molecule = new " + molecule + " [1]\n\n")
    system.write('write_once("Data Boundary") {\n')
    system.write("\t0.0 " + box + " xlo xhi\n")
    system.write("\t0.0 " + box + " ylo yhi\n")
    system.write("\t0.0 " + box + " zlo zhi\n")
    system.write("}\n\n")
