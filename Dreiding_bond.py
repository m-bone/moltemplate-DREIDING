# Import packages
import os
import pandas as pd
import numpy as np
from Dreiding_label_dictionary import labelDict

os.chdir("/home/matt/Documents/Dreiding_forcefield")

subDict = labelDict

# Replace many labels with a single wildcard label
# Wildcard function drops values and replaces them with replacement
def wildcard(list, replacement):
    for value in list:
        subDict[replacement] = subDict.pop(value)

# Remove duplicate labels that don't effect bonding
# Carbon
wildcard(["C_3", "C_34", "C_33", "C_32", "C_31"], "C_3*")
wildcard(["C_R", "C_R1", ], "C_R*")
wildcard(["C_R_b1", "C_R1_b1"], "C_R*_b1")

# Nitrogen
wildcard(["N_3", "N_3_ha", "N_3_hd"], "N_3*")
wildcard(["N_R_d1", "N_R_d1_ha", "N_R_d2", "N_R_d2_ha", "N_R_d2_hd"], "N_R*")
wildcard(["N_R_b1_d2", "N_R_b1_d2_ha", "N_R_b1_d2_hd"], "N_R_b1*")
wildcard(["N_2_d1", "N_2_d1_ha", "N_2_d1_hd", "N_2_d2", "N_2_d2_ha", "N_2_d2_hd"], "N_2*")
wildcard(["N_2_b1_d1", "N_2_b1_d1_ha","N_2_b1_d2", "N_2_b1_d2_ha", "N_2_b1_d2_hd"], "N_2_b1*")
wildcard(["N_2_b2_d1", "N_2_b2_d1_ha", "N_2_b2_d2", "N_2_b2_d2_ha", "N_2_b2_d2_hd"], "N_2_b2*")
wildcard(["N_1", "N_1_ha"], "N_1*")

# Oxygen
wildcard(["O_3", "O_3_ha", "O_3_hd"], "O_3*")
wildcard(["O_R", "O_R_ha"], "O_R*")
wildcard(["O_2", "O_2_ha"], "O_2*")
wildcard(["O_2_b1", "O_2_b1_ha", "O_2_b1_hd"], "O_2_b1*")
wildcard(["O_2_b2", "O_2_b2_ha"], "O_2_b2*")
wildcard(["O_1", "O_1_ha"], "O_1*")

# Boron
wildcard(["B_2_d1", "B_2_d2", ], "B_2*")
wildcard(["B_2_b1_d1", "B_2_b1_d2", ], "B_2_b1*")
wildcard(["B_2_b2_d1", "B_2_b2_d2", ], "B_2_b2*")

# Fluorine
wildcard(["F", "F_ha", "F_hd"], "F*")

# Other elements
wildcard(["Al_3_d1", "Al_3_d2"], "Al_3*")
wildcard(["P_3_d2", "P_3_d3", "P_3_d4"], "P_3*")
wildcard(["Ge_3_d1", "Ge_3_d2", "Ge_3_d3"], "Ge_3*")
wildcard(["As_3_d1", "As_3_d2", "As_3_d3", "As_3_d4"], "As_3*")
wildcard(["Se_3_d1", "Se_3_d2", "Se_3_d3", "Se_3_d5"], "Se_3*")
wildcard(["In_3_d1", "In_3_d2"], "In_3*")
wildcard(["Sn_3_d2", "Sn_3_d3"], "Sn_3*")
wildcard(["Sb_3_d2", "Sb_3_d3", "Sb_3_d4"], "Sb_3*")
wildcard(["Te_3_d1", "Te_3_d2", "Te_3_d3", "Te_3_d5"], "Te_3*")

#List of keys from dictionary
keys = list(subDict.keys())

# Dictionary of bond radii according to DREIDING
atomClasses = {key: subDict[key][1] for key in subDict.keys()}

# Empty dataframe with rows and columns labelled as keys
R0 = pd.DataFrame([], index = [keys], columns = [keys])


# R0 calculation from DREIDING
for i in range(len(atomClasses)):
    for j in range(len(atomClasses)):
        R0.loc[keys[i],keys[j]] = atomClasses[keys[i]] + atomClasses[keys[j]] - 0.01


# Convert pandas dataframe to numpy matrix so the diagonal can be found
R0 = R0.to_numpy()
R0 = R0.astype('float64')
R0 = np.asmatrix(R0)
R0 = np.tril(R0)
R0 = np.around(R0, decimals = 3)

# For debug checking
R0Pandas = pd.DataFrame(R0)
R0Pandas.columns = keys
R0Pandas.index = keys
R0Pandas

# Set K in (kcal/mol)/A**2
K = 700/2

# Create files for bond by type and bond coeff
file_type = open("bondsType.txt", "w")
file_coeff = open("bondsCoeff.txt", "w")

# Store of previously formed bondNames
def bondWrite(multiplier):
        bondName = keys[i] + "-" + keys[j]
        bondName = bondName.replace("*", "w")

        # Write to files and append bondName to prevBond
        file_type.write("@bond:" + bondName + " @atom:" + keys[i] + " @atom:" + keys[j] + "\n")
        file_coeff.write("bond_coeff @bond:" + bondName + " " + str(K*multiplier) + " " + str(R0[i, j]) + "\n")


for i in range(len(atomClasses)):
    for j in range(len(atomClasses)):
        if R0[i,j] == 0: # Skip zero value cells in matrix
            continue

        # Catch X_b1-x_b1 bonds
        if "_b1" in keys[i] and "_b1" in keys[j]:
            bondWrite(1)

        # Catch X_2_b2-X_2_b2 bonds where bond order is 2
        elif "_b2" in keys[i] and "_b2" in keys[j]:
            bondWrite(2)

        # Catch X_2_b2-X_2 bonds where bond order is 1
        elif "_b2" in keys[i] and "_2" in keys[j]:
            bondWrite(1)
        elif "_2" in keys[i] and "_b2" in keys[j]:
            bondWrite(1)

        # Catch X_1-X_1 bonds where the bond order will always be 3
        elif "_1" in keys[i] and "_1" in keys[j]:
            bondWrite(3)

        # Catch X_2-X_2 bonds where the bond order will always be 2
        elif "_2" in keys[i] and "_2" in keys[j]:
            bondWrite(2)

        # Catch X_R-X_R bonds where the bond order will always be 1.5
        elif "_R" in keys[i] and "_R" in keys[j]:
            bondWrite(1.5)

        # If not a specialist case, handle writing in a general fashion
        else:
            bondWrite(1)

file_type.close()
file_coeff.close()
