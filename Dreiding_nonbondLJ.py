##############################################################################
# Developed by: Matthew Bone
# Last Updated: 01/10/2020
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
# This file generates moltemplate coefficient data for the nonbonding Lennard -
# Jones potential. As DREIDING only uses element to determine nonbonding
# coefficients, this script uses an atomDictionary rather than the
# Dreiding_label_dictionary. All R values are the original DREIDING value
# divided by 2^1/6 as LAMMPS requires R at the zero crossing point.
##############################################################################
# Import packages
import os
import pandas as pd
import numpy as np

os.chdir("/home/matt/Documents/Dreiding_forcefield")

# [0] is R, [1] is D
atomDictionary = {
    "H": [2.846421, 0.0152],
    "H_HB": [2.846421, 0.0001],
    "H_B": [2.846421, 0.0152],
    "B*": [3.581513, 0.095],
    "C*": [3.472990, 0.0951],
    "N*": [3.262560, 0.0774],
    "O*": [3.033153, 0.0957],
    "F*": [3.093200, 0.0725],
    "Al*": [3.893227, 0.31],
    "Si*": [3.804138, 0.31],
    "P*": [3.697230, 0.32],
    "S*": [3.590322, 0.32],
    "Cl": [3.519317, 0.2833],
    "Ga*": [3.911045, 0.40],
    "Ge*": [3.804138, 0.40],
    "As*": [3.697230, 0.41],
    "Se*": [3.590321, 0.43],
    "Br": [3.519050, 0.37],
    "In*": [4.089225, 0.55],
    "Sn*": [3.982317, 0.55],
    "Sb*": [3.875409, 0.55],
    "Te*": [3.768502, 0.57],
    "I": [3.697230, 0.51],
    "Na": [2.800986, 0.50],
    "Ca": [3.093200, 0.050],
    "Fe": [4.044680, 0.055],
    "Zn": [4.044680, 0.055],

    "C_R1": [3.768502, 0.1356],
    "C_34": [3.774738, 0.3016],
    "C_33": [3.699368, 0.2500],
    "C_32": [3.623909, 0.1984],
    "C_31": [3.548450, 0.1467],
}

# List of keys from dictionary
keys = list(atomDictionary.keys())

# Empty dataframe for Rij with rows and columns labelled as keys
Rij = pd.DataFrame([], index = [keys], columns = [keys])

# Rij calculation from DREIDING
for i in range(len(atomDictionary)):
    for j in range(len(atomDictionary)):
        Rij.loc[keys[i],keys[j]] = (atomDictionary[keys[i]][0] + atomDictionary[keys[j]][0]) / 2

# Convert to numpy matrix so that values below the diagonal can be found
Rij = Rij.to_numpy()
Rij = Rij.astype('float64')
Rij = np.asmatrix(Rij)
Rij = np.tril(Rij)
Rij = np.around(Rij, decimals = 8)
Rij

# Empty dataframe for Dij with rows and columns labelled as keys
Dij = pd.DataFrame([], index = [keys], columns = [keys])

# Dij calculation from DREIDING
for i in range(len(atomDictionary)):
    for j in range(len(atomDictionary)):
        Dij.loc[keys[i],keys[j]] = (atomDictionary[keys[i]][1] * atomDictionary[keys[j]][1])**0.5

# Convert to numpy matrix so that values below the diagonal can be found
Dij = Dij.to_numpy()
Dij = Dij.astype('float64')
Dij = np.asmatrix(Dij)
Dij = np.tril(Dij)
Dij = np.around(Dij, decimals = 8)
Dij

# Create text file
file = open("lj_nonbond.txt", "w")

# Write lines of text file using wildcards for C, N and O; not H and H_HB
for i in range(len(atomDictionary)):
    for j in range(len(atomDictionary)):
        if Rij[i,j] == 0:
            continue
        file.write("pair_coeff  @atom:" + keys[i] + "  @atom:" + keys[j] + " lj/cut/coul/long " + str(Dij[i,j]) + " " + str(Rij[i,j]) + "\n")

file.close()
