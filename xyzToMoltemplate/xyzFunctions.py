##############################################################################
# Developed by: Matthew Bone
# Last Updated: 08/10/2020
# Updated by: Matthew Bone
#
# Contact Details:
# Bristol Composites Institute (BCI)
# Department of Aerospace Engineering - University of Bristol
# Queen's Building - University Walk
# Bristol, BS8 1TR
# U.K.
# Email - matthew.bone@bristol.ac.uk
#
# File Description:
# This file contains functions to manipulate .xyz files and produce moltemplate
# .lt files. It is typically used with xyzToMoltemplate.py. The atom dictionary
# can be updated to included more atoms and masses.
##############################################################################

# Import packages
import math
import numpy as np
import pandas as pd

# Dictionary of relative atomic mass
# Values taken from https://www.rsc.org/periodic-table
atomDictionary = {
    "H": 1.008,
    "He": 4.003,
    "Li": 6.940,
    "Be": 9.012,
    "B": 10.810,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "Ne": 20.180,
    "Na": 22.990,
    "Mg": 24.305,
    "Al": 26.982,
    "Si": 28.085,
    "P": 30.974,
    "S": 32.060,
    "Cl": 35.450,
    "Ar": 39.950,
    #...
    "Ca": 40.078,
    "Fe": 55.845,
    "Zn": 65.380,
    "Ga": 69.723,
    "Ge": 72.630,
    "As": 74.922,
    "Se": 78.971,
    "Br": 79.904,
    "In": 114.818,
    "Sn": 118.710,
    "Sb": 121.760,
    "Te": 127.600,
    "I": 126.904,
}

# Import xyz file and store as pd.DataFrame
def xyzToDF(filename):
    # Import xyz file using np
    atomData = np.genfromtxt(filename, delimiter=None, dtype=[("Element", '<U32'),
     ("X", float), ("Y", float), ("Z", float)], skip_header=2, encoding="UTF8")
    # Convert to pd dataframe
    atomData = pd.DataFrame(atomData)

    return atomData

# Determine what atoms are bonded by looking at bonding distance
def bondCalculation(atomData):
    # Preload list for storing bond data
    bondData = []

    # Create a list of tuples of (int, list) that contain the bonder (bondData[x][0]) and bondees (bondData[x][1][y])
    for index1, row1 in atomData.iterrows():
        if row1["Element"] == 'H':
            continue # Prevents creating bonds with H as the bonder. Not needed as H can only ever have one bond
        bondList = []
        for index2, row2 in atomData.iterrows():
            if index1 == index2:
                continue # Stop bonder forming bonds with itself
            bondLength = math.sqrt((row1["X"] - row2["X"])**2 + (row1["Y"] - row2["Y"])**2 + (row1["Z"] - row2["Z"])**2) # Calc bond length
            if bondLength <= 1.6: # As structures are very ordered. Anything within 1.6 Ang is going to be a bonded atom
                bondList.append(index2)
        bondData.append((index1, bondList))

    return bondData


# Add moletemplate jargon to column including element and element count
def moltemplater(atomData, atomDictionary):
    # Preload list for moltemplate atom names
    atomName = []

    # Dataframe of elements and counters
    countersDF = pd.DataFrame({"Counter": 0}, index = list(atomDictionary.keys()))

    # Check to see if there are elements in the file that don't match the database
    incomingElements = set(list(atomData["Element"]))
    databaseElements = set(list(atomDictionary.keys()))
    differentElements = incomingElements.difference(databaseElements)

    # Error and return unmatched elements
    assert(len(incomingElements.intersection(databaseElements)) == len(incomingElements)), "File contains unknown elements: {}".format(differentElements)

    # Loop through atomData
    for index, row in atomData.iterrows():
        # Find element symbol and add moltemplate syntax, element symbol and counter value to atomData
        atom = row["Element"]
        atomData.loc[index, "Element"] = "$atom:" + row["Element"] + str(countersDF.loc[atom, "Counter"])

        # Add atom name to atomName list
        atomName.append(atom + str(countersDF.loc[atom, "Counter"]))

        # Increase relevant counter by 1
        countersDF.loc[atom] += 1

    return atomName
