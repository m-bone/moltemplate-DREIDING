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
# This file generates moltemplate coefficient data for the dihedral angle section
# of the force field. Labels are taken from the Dreiding_label_dictionary. The
# dihedral calculation is the most complicated of all the force field components
# as it relies on 10 different cases to catagorise all possible dihedral
# coefficients. As such this file may take around 30 seconds to compute.
#
# All cases are built by subsetting the main label dictionary.
# The wildcard functions compress labels with unnecessary flags into labels with
# a wildcard (*). This is done later in some case dictionaries as flags are
# required to subset the main dictionary. Constants can be changed for each case
# at the start of each case. All DREIDING dihedral cases only use the dihedral
# bonding atoms, apart from Case J. No dihedral type names should be repeated;
# there is a check to generate an error if this is happens. Case I is at the
# bottom so that it has priority when computed by Moltemplate
##############################################################################

# Import packages
import os
from Dreiding_label_dictionary import labelDict

os.chdir("/home/matt/Documents/Dreiding_forcefield")
subDict = labelDict

# List of keys from dictionary
keys = list(subDict.keys())

# Replace many labels with a single wildcard label
# Wildcard function drops values and replaces them with replacement
def wildcard(list, replacement):
    for value in list:
        subDict[replacement] = subDict.pop(value)

def wildcardCases(list, replacement, dict):
    for value in list:
        dict[replacement] = dict.pop(value)

# Reducing labels with replacement wildcard - wildcard([], )
# Hydrogen
wildcard(["H", "H_HB", "H_B"], "H*")

# Carbon
wildcard(["C_3", "C_34", "C_33", "C_32", "C_31"], "C_3*")
wildcard(["C_R", "C_R1", ], "C_R*")
wildcard(["C_R_b1", "C_R1_b1"], "C_R*_b1")
# wildcard(["C_2", "C_2_b1", "C_2_b2"], "C_2*")
wildcard(["C_1", "C_1_b1"], "C_1*")

# Nitrogen
wildcard(["N_3", "N_3_ha", "N_3_hd"], "N_3*")
wildcard(["N_R_d1", "N_R_d1_ha"], "N_R_d1*")
wildcard(["N_R_d2", "N_R_d2_ha", "N_R_d2_hd"], "N_R_d2*")
wildcard(["N_R_b1_d2", "N_R_b1_d2_ha", "N_R_b1_d2_hd"], "N_R_b1_d2*")
wildcard(["N_2_d1", "N_2_d1_ha", "N_2_d1_hd"], "N_2_d1*")
wildcard(["N_2_b1_d1", "N_2_b1_d1_ha"], "N_2_b1_d1*")
wildcard(["N_2_b2_d1", "N_2_b2_d1_ha"], "N_2_b2_d1*")
wildcard(["N_2_d2", "N_2_d2_ha", "N_2_d2_hd"], "N_2_d2*")
wildcard(["N_2_b1_d2", "N_2_b1_d2_ha", "N_2_b1_d2_hd"], "N_2_b1_d2*")
wildcard(["N_2_b2_d2", "N_2_b2_d2_ha", "N_2_b2_d2_hd"], "N_2_b2_d2*")
wildcard(["N_1", "N_1_ha"], "N_1*")

# Oxygen
wildcard(["O_3", "O_3_ha", "O_3_hd"], "O_3*")
wildcard(["O_R", "O_R_ha"], "O_R*")
wildcard(["O_2", "O_2_ha"], "O_2*")
wildcard(["O_2_b1", "O_2_b1_ha", "O_2_b1_hd"], "O_2_b1*")
wildcard(["O_2_b2", "O_2_b2_ha"], "O_2_b2*")
wildcard(["O_1", "O_1_ha"], "O_1*")


# Fluorine
wildcard(["F", "F_hd", "F_ha"], "F*")

# Number of possible atomic bonds, aside from the J-K bond in question
possibleAtoms = {key: subDict[key][3] for key in subDict.keys()}

# List of X_3 in column 16 as treated differently (oxygen, sulphur, etc.)
column16 = ["O_3*", "S_3", "Se_3", "Te_3"]

# Repeated type names list: No type names should appear twice with this methodology
typeNames = []

# Start dihedral by type and coeff files
file_type = open("dihedralType.txt", "w")
file_coeff = open("dihedralCoeff.txt", "w")

# Split into dihedral cases from DREIDING paper
# Case A - J,K = X_3
# Subpriority with Case H
CASE_A_ENERGY = 2.0/2
CASE_A_MULTIPLICITY = 3
CASE_A_PHASE_SHIFT = 180 * CASE_A_MULTIPLICITY + 180

# Subsets main dictionary for values with "_3" in the key
caseADictionary = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_3" in key}
caseAKeys = list(caseADictionary.keys())

for j in range(len(caseADictionary)):
    for k in range(len(caseADictionary)):

        # Due to Case H priority - skips column16-column16 bonding pairs
        if caseAKeys[j][0:4] in column16 and caseAKeys[k][0:4] in column16:
            continue

        # Create type name for dihedral, remove wildcard * if present
        dihedralName = caseAKeys[j] + "-" + caseAKeys[k]
        dihedralName = dihedralName.replace("*", "w")

        # Writes to files, calculate barrier energy and add key combo to list
        file_type.write("@dihedral:" + dihedralName + "  @atom:* @atom:" + caseAKeys[j] + " @atom:" + caseAKeys[k] + " @atom:*\n")
        barrierEnergy = str(round(CASE_A_ENERGY / (caseADictionary[caseAKeys[j]] * caseADictionary[caseAKeys[k]]), 4))
        file_coeff.write("dihedral_coeff @dihedral:" + dihedralName + " " + barrierEnergy + " " + str(CASE_A_MULTIPLICITY) + " " + str(CASE_A_PHASE_SHIFT) + " 0.000\n")
        typeNames.append(dihedralName)

# Case B - J = X_2, X_R, K = X_3 single bond with one sp2 and sp3 atoms e.g. toluene
CASE_B_ENERGY = 1.0/2
CASE_B_MULTIPLICITY = 6
CASE_B_PHASE_SHIFT = 0 * CASE_B_MULTIPLICITY + 180

# Subsets main dictionary for values with "_3", "_2" or "_R" in the key
caseBDictionaryJ = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_2" in key or "_R" in key}
wildcardCases(["C_2", "C_2_b1", "C_2_b2"], "C_2*", caseBDictionaryJ)
wildcardCases(["C_R*", "C_R*_b1"], "C_R*", caseBDictionaryJ)
wildcardCases(["B_2_d1", "B_2_b1_d1", "B_2_b2_d1"], "B_2*_d1", caseBDictionaryJ)
wildcardCases(["B_2_d2", "B_2_b1_d2", "B_2_b2_d2"], "B_2*_d2", caseBDictionaryJ)
wildcardCases(["N_R_d2*", "N_R_b1_d2*"], "N_R*_d2*", caseBDictionaryJ)
wildcardCases(["N_2_d1*", "N_2_b1_d1*", "N_2_b2_d1*"], "N_2*_d1*", caseBDictionaryJ)
wildcardCases(["N_2_d2*", "N_2_b1_d2*", "N_2_b2_d2*"], "N_2*_d2*", caseBDictionaryJ)
wildcardCases(["O_2*", "O_2_b1*", "O_2_b2*"], "O_2*", caseBDictionaryJ)

caseBDictionaryK = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_3" in key}
caseBKeysK = list(caseBDictionaryK.keys())
caseBKeysJ = list(caseBDictionaryJ.keys())

for j in range(len(caseBDictionaryJ)):
    for k in range(len(caseBDictionaryK)):

        # Due to Case I priority - Case J given priority by MT rather than here
        if keys[k] in column16:
            continue

        # Create type name for dihedral, remove wildcard * if present
        dihedralName = caseBKeysJ[j] + "-" + caseBKeysK[k]
        dihedralName = dihedralName.replace("*", "w")

        # Writes to files, calculate barrier energy and add key combo to list
        file_type.write("@dihedral:" + dihedralName + "  @atom:* @atom:" + caseBKeysJ[j] + " @atom:" + caseBKeysK[k] + " @atom:*\n")
        barrierEnergy = str(round(CASE_B_ENERGY / (caseBDictionaryJ[caseBKeysJ[j]] * caseBDictionaryK[caseBKeysK[k]]), 4))
        file_coeff.write("dihedral_coeff @dihedral:" + dihedralName + " " + barrierEnergy + " " + str(CASE_B_MULTIPLICITY) + " " + str(CASE_B_PHASE_SHIFT) + " 0.000\n")
        typeNames.append(dihedralName)

# Case C - J,K = X_2 double bond with two sp2 atoms
CASE_C_ENERGY = 45.0/2
CASE_C_MULTIPLICITY = 2
CASE_C_PHASE_SHIFT = 180 * CASE_C_MULTIPLICITY + 180

# Subsets main dictionary for values with "_2" in the key
# Two dictionaries as need one with _b1 and one without
caseCDictionaryJ = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_2" in key}
caseCDictionaryK = {key: caseCDictionaryJ[key] for key in caseCDictionaryJ.keys() if "_b1" not in key}
caseCKeysJ = list(caseCDictionaryJ.keys())
caseCKeysK = list(caseCDictionaryK.keys())


for j in range(len(caseCDictionaryJ)):
    for k in range(len(caseCDictionaryK)):

        # Removes "notb2"-b2 bonds which are single bonds
        if "_b2" in caseCKeysJ[j] and "_b2" not in caseCKeysK[k]:
            continue
        if "_b2" not in caseCKeysJ[j] and "_b2" in caseCKeysK[k]:
            continue

        # Create type name for dihedral, remove wildcard * if present
        dihedralName = caseCKeysJ[j] + "-" + caseCKeysK[k]
        dihedralName = dihedralName.replace("*", "w")

        # Writes to files, calculate barrier energy and add key combo to list
        file_type.write("@dihedral:" + dihedralName + "  @atom:* @atom:" + caseCKeysJ[j] + " @atom:" + caseCKeysK[k] + " @atom:*\n")
        barrierEnergy = str(round(CASE_C_ENERGY / (caseCDictionaryJ[caseCKeysJ[j]] * caseCDictionaryK[caseCKeysK[k]]), 4))
        file_coeff.write("dihedral_coeff @dihedral:" + dihedralName + " " + barrierEnergy + " " + str(CASE_C_MULTIPLICITY) + " " + str(CASE_C_PHASE_SHIFT) + " 0.000\n")
        typeNames.append(dihedralName)

# Case D - J,K = X_R resonance bond involving two resonance atoms
CASE_D_ENERGY = 25.0/2
CASE_D_MULTIPLICITY = 2
CASE_D_PHASE_SHIFT = 180 * CASE_D_MULTIPLICITY + 180

# Subsets main dictionary for values with "_R" in the key
# Two dictionaries as need one with _b1 and one without
caseDDictionaryJ = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_R" in key}
caseDDictionaryK = {key: caseDDictionaryJ[key] for key in caseDDictionaryJ.keys() if "_b1" not in key}
caseDKeysJ = list(caseDDictionaryJ.keys())
caseDKeysK = list(caseDDictionaryK.keys())

for j in range(len(caseDDictionaryJ)):
    for k in range(len(caseDDictionaryK)):

        # Create type name for dihedral, remove wildcard * if present
        dihedralName = caseDKeysJ[j] + "-" + caseDKeysK[k]
        dihedralName = dihedralName.replace("*", "w")

        # Writes to files, calculate barrier energy and add key combo to list
        file_type.write("@dihedral:" + dihedralName + "  @atom:* @atom:" + caseDKeysJ[j] + " @atom:" + caseDKeysK[k] + " @atom:*\n")
        barrierEnergy = str(round(CASE_D_ENERGY / (caseDDictionaryJ[caseDKeysJ[j]] * caseDDictionaryK[caseDKeysK[k]]), 4))
        file_coeff.write("dihedral_coeff @dihedral:" + dihedralName + " " + barrierEnergy + " " + str(CASE_D_MULTIPLICITY) + " " + str(CASE_D_PHASE_SHIFT) + " 0.000\n")
        typeNames.append(dihedralName)

# Case E - J,K = X_2, X_R single bond involving two sp2 or resonant atoms ONLY B1
CASE_E_ENERGY = 5.0/2
CASE_E_MULTIPLICITY = 2
CASE_E_PHASE_SHIFT = 180 * CASE_E_MULTIPLICITY + 180

# Subsets main dictionary for values with "_R_b1" or "_2_b1" in the key
# Two dictionaries so that R_b1-R_b1 bonds are left for Case F
# Split into two halves due to complexity added with the b2 flag
caseEDictionaryJ = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_R_b1" in key or "_R*_b1" in key or "_2_b1" in key}
caseEDictionaryK = {key: caseEDictionaryJ[key] for key in caseEDictionaryJ.keys() if "_R_b1" not in key}
caseEDictionaryK.pop("C_R*_b1")
caseEKeysJ = list(caseEDictionaryJ.keys())
caseEKeysK = list(caseEDictionaryK.keys())


for j in range(len(caseEDictionaryJ)):
    for k in range(len(caseEDictionaryK)):

        # Create type name for dihedral, remove wildcard * if present
        dihedralName = caseEKeysJ[j] + "-" + caseEKeysK[k]
        dihedralName = dihedralName.replace("*", "w")

        # Writes to files, calculate barrier energy and add key combo to list
        file_type.write("@dihedral:" + dihedralName + "  @atom:* @atom:" + caseEKeysJ[j] + " @atom:" + caseEKeysK[k] + " @atom:*\n")
        barrierEnergy = str(round(CASE_E_ENERGY / (caseEDictionaryJ[caseEKeysJ[j]] * caseEDictionaryK[caseEKeysK[k]]), 4))
        file_coeff.write("dihedral_coeff @dihedral:" + dihedralName + " " + barrierEnergy + " " + str(CASE_E_MULTIPLICITY) + " " + str(CASE_E_PHASE_SHIFT) + " 0.000\n")
        typeNames.append(dihedralName)

# Second half for catching b2-other bonds
caseEDictionaryJ = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "b2" in key}
caseEDictionaryK = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_2" in key and "b2" not in key}
caseEKeysJ = list(caseEDictionaryJ.keys())
caseEKeysK = list(caseEDictionaryK.keys())

for j in range(len(caseEDictionaryJ)):
    for k in range(len(caseEDictionaryK)):

        # Create type name for dihedral, remove wildcard * if present
        dihedralName = caseEKeysJ[j] + "-" + caseEKeysK[k]
        dihedralName = dihedralName.replace("*", "w")

        # Writes to files, calculate barrier energy and add key combo to list
        file_type.write("@dihedral:" + dihedralName + "  @atom:* @atom:" + caseEKeysJ[j] + " @atom:" + caseEKeysK[k] + " @atom:*\n")
        barrierEnergy = str(round(CASE_E_ENERGY / (caseEDictionaryJ[caseEKeysJ[j]] * caseEDictionaryK[caseEKeysK[k]]), 4))
        file_coeff.write("dihedral_coeff @dihedral:" + dihedralName + " " + barrierEnergy + " " + str(CASE_E_MULTIPLICITY) + " " + str(CASE_E_PHASE_SHIFT) + " 0.000\n")
        typeNames.append(dihedralName)

# Second half for catching b2-other bonds
caseEDictionaryJ = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_R" in key and "_b1" not in key}
caseEDictionaryK = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_2" in key}
caseEKeysJ = list(caseEDictionaryJ.keys())
caseEKeysK = list(caseEDictionaryK.keys())

for j in range(len(caseEDictionaryJ)):
    for k in range(len(caseEDictionaryK)):

        # Create type name for dihedral, remove wildcard * if present
        dihedralName = caseEKeysJ[j] + "-" + caseEKeysK[k]
        dihedralName = dihedralName.replace("*", "w")

        # Writes to files, calculate barrier energy and add key combo to list
        file_type.write("@dihedral:" + dihedralName + "  @atom:* @atom:" + caseEKeysJ[j] + " @atom:" + caseEKeysK[k] + " @atom:*\n")
        barrierEnergy = str(round(CASE_E_ENERGY / (caseEDictionaryJ[caseEKeysJ[j]] * caseEDictionaryK[caseEKeysK[k]]), 4))
        file_coeff.write("dihedral_coeff @dihedral:" + dihedralName + " " + barrierEnergy + " " + str(CASE_E_MULTIPLICITY) + " " + str(CASE_E_PHASE_SHIFT) + " 0.000\n")
        typeNames.append(dihedralName)

# Case F - J,K = X_R
# Lower so greater priority over case E, which it supersedes in some cases
CASE_F_ENERGY = 10.0/2
CASE_F_MULTIPLICITY = 2
CASE_F_PHASE_SHIFT = 180 * CASE_F_MULTIPLICITY + 180

# Subsets main dictionary for values with "_R_b1" in the key
# This handles R_b1-R_b1 bonds
caseFDictionary = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_R_b1" in key or "_R*_b1" in key}
caseFKeys = list(caseFDictionary.keys())

for j in range(len(caseFDictionary)):
    for k in range(len(caseFDictionary)):

        # Create type name for dihedral, remove wildcard * if present
        dihedralName = caseFKeys[j] + "-" + caseFKeys[k]
        dihedralName = dihedralName.replace("*", "w")

        # Writes to files, calculate barrier energy and add key combo to list
        file_type.write("@dihedral:" + dihedralName + "  @atom:* @atom:" + caseFKeys[j] + " @atom:" + caseFKeys[k] + " @atom:*\n")
        barrierEnergy = str(round(CASE_F_ENERGY / (caseFDictionary[caseFKeys[j]] * caseFDictionary[caseFKeys[k]]), 4))
        file_coeff.write("dihedral_coeff @dihedral:" + dihedralName + " " + barrierEnergy + " " + str(CASE_F_MULTIPLICITY) + " " + str(CASE_F_PHASE_SHIFT) + " 0.000\n")
        typeNames.append(dihedralName)

# Case G - J,K = X_1, monovalent, metal
# This case is included to make up the numbers in terms of dihderals
# Uses wildcards to prevent significant addition of lines to force field
CASE_G_ENERGY = 0
CASE_G_MULTIPLICITY = 0
CASE_G_PHASE_SHIFT = 0

# Subsets main dictionary for values with "_1" in the key
caseGDictionary = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_1" in key}
caseGKeys = list(caseGDictionary.keys())

for j in range(len(caseGDictionary)):

    # Writes to files, calculate barrier energy and add key combo to list
    file_type.write("@dihedral:" + caseGKeys[j].replace("*", "w") + "-wildcard" + "  @atom:* @atom:" + caseGKeys[j] + " @atom:*" + " @atom:*\n")
    file_coeff.write("dihedral_coeff @dihedral:" + caseGKeys[j].replace("*", "w") + "-wildcard " + str(CASE_G_ENERGY) + " " + str(CASE_G_MULTIPLICITY) + " " + str(CASE_G_PHASE_SHIFT) + " 0.000\n")
    typeNames.append(caseGKeys[j].replace("*", "w") + "-wildcard")

# Case H - J,K = X_3 from column 16
# Alternate form of Case A
CASE_H_ENERGY = 2.0/2
CASE_H_MULTIPLICITY = 2
CASE_H_PHASE_SHIFT = 90 * CASE_H_MULTIPLICITY + 180

# Subsets main dictionary for values with "_3" in the key if they are column16
caseHDictionary = {key: possibleAtoms[key] for key in possibleAtoms.keys() for x3 in column16 if x3 in key}
caseHKeys = list(caseHDictionary.keys())

for j in range(len(caseHDictionary)):
    for k in range(len(caseHDictionary)):

        # Create type name for dihedral, remove wildcard * if present
        dihedralName = caseHKeys[j] + "-" + caseHKeys[k]
        dihedralName = dihedralName.replace("*", "w")

        # Writes to files, calculate barrier energy and add key combo to list
        file_type.write("@dihedral:" + dihedralName + "  @atom:* @atom:" + caseHKeys[j] + " @atom:" + caseHKeys[k] + " @atom:*\n")
        barrierEnergy = str(round(CASE_H_ENERGY / (caseHDictionary[caseHKeys[j]] * caseHDictionary[caseHKeys[k]]), 4))
        file_coeff.write("dihedral_coeff @dihedral:" + dihedralName + " " + barrierEnergy + " " + str(CASE_H_MULTIPLICITY) + " " + str(CASE_H_PHASE_SHIFT) + " 0.000\n")
        typeNames.append(dihedralName)

# Case I - J = X_3 from column 16, K = X_2, X_R
# Alternate form of Case B
CASE_I_ENERGY = 2.0/2
CASE_I_MULTIPLICITY = 2
CASE_I_PHASE_SHIFT = 180 * CASE_I_MULTIPLICITY + 180

# Subsets main dictionary for values with "_3" for column16, "_2" or "_R" in the key
# C_R_b1 can't techically ever do this but has been left
caseIDictionaryJ = {key: possibleAtoms[key] for key in possibleAtoms.keys() for x3 in column16 if x3 in key}
caseIDictionaryK = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_2" in key or "_R" in key}
wildcardCases(["C_2", "C_2_b1", "C_2_b2"], "C_2*", caseIDictionaryK)
wildcardCases(["C_R*", "C_R*_b1"], "C_R*", caseIDictionaryK)
wildcardCases(["B_2_d1", "B_2_b1_d1", "B_2_b2_d1"], "B_2*_d1", caseIDictionaryK)
wildcardCases(["B_2_d2", "B_2_b1_d2", "B_2_b2_d2"], "B_2*_d2", caseIDictionaryK)
wildcardCases(["N_R_d2*", "N_R_b1_d2*"], "N_R*_d2*", caseIDictionaryK)
wildcardCases(["N_2_d1*", "N_2_b1_d1*", "N_2_b2_d1*"], "N_2*_d1*", caseIDictionaryK)
wildcardCases(["N_2_d2*", "N_2_b1_d2*", "N_2_b2_d2*"], "N_2*_d2*", caseIDictionaryK)
wildcardCases(["O_2*", "O_2_b1*", "O_2_b2*"], "O_2*", caseIDictionaryK)

caseIKeysK = list(caseIDictionaryK.keys())
caseIKeysJ = list(caseIDictionaryJ.keys())

for j in range(len(caseIDictionaryJ)):
    for k in range(len(caseIDictionaryK)):

        # Create type name for dihedral, remove wildcard * if present
        dihedralName = caseIKeysJ[j] + "-" + caseIKeysK[k]
        dihedralName = dihedralName.replace("*", "w")

        # Writes to files, calculate barrier energy and add key combo to list
        file_type.write("@dihedral:" + dihedralName + "  @atom:* @atom:" + caseIKeysJ[j] + " @atom:" + caseIKeysK[k] + " @atom:*\n")
        barrierEnergy = str(round(CASE_I_ENERGY / (caseIDictionaryJ[caseIKeysJ[j]] * caseIDictionaryK[caseIKeysK[k]]), 4))
        file_coeff.write("dihedral_coeff @dihedral:" + dihedralName + " " + barrierEnergy + " " + str(CASE_I_MULTIPLICITY) + " " + str(CASE_I_PHASE_SHIFT) + " 0.000\n")
        typeNames.append(dihedralName)

# Case J - J= X_2, X_R, K = X_3, I != X_2, X_R e.g. propene
# Alternate form of Case B - Last means it has greatest priority in moltemplate
CASE_J_ENERGY = 2.0/2
CASE_J_MULTIPLICITY = 3
CASE_J_PHASE_SHIFT = 180 * CASE_J_MULTIPLICITY + 180

caseJDictionaryJ = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_2" in key or "_R" in key}
wildcardCases(["C_2", "C_2_b1", "C_2_b2"], "C_2*", caseJDictionaryJ)
wildcardCases(["C_R*", "C_R*_b1"], "C_R*", caseJDictionaryJ)
wildcardCases(["B_2_d1", "B_2_b1_d1", "B_2_b2_d1"], "B_2*_d1", caseJDictionaryJ)
wildcardCases(["B_2_d2", "B_2_b1_d2", "B_2_b2_d2"], "B_2*_d2", caseJDictionaryJ)
wildcardCases(["N_R_d2*", "N_R_b1_d2*"], "N_R*_d2*", caseJDictionaryJ)
wildcardCases(["N_2_d1*", "N_2_b1_d1*", "N_2_b2_d1*"], "N_2*_d1*", caseJDictionaryJ)
wildcardCases(["N_2_d2*", "N_2_b1_d2*", "N_2_b2_d2*"], "N_2*_d2*", caseJDictionaryJ)
wildcardCases(["O_2*", "O_2_b1*", "O_2_b2*"], "O_2*", caseJDictionaryJ)

caseJDictionaryK = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_3" in key}
caseJDictionaryI = {key: possibleAtoms[key] for key in possibleAtoms.keys() if "_2" not in key and "_R" not in key}
caseJKeysK = list(caseJDictionaryK.keys())
caseJKeysJ = list(caseJDictionaryJ.keys())
caseJKeysI = list(caseJDictionaryI.keys())

for j in range(len(caseJDictionaryJ)):
    for k in range(len(caseJDictionaryK)):
        for i in range(len(caseJDictionaryI)):

            # Create type name for dihedral, remove wildcard * if present
            dihedralName = caseJKeysJ[j] + "-" + caseJKeysK[k] + "-" + caseJKeysI[i]
            dihedralName = dihedralName.replace("*", "w")

            # Writes to files, calculate barrier energy and add key combo to list
            file_type.write("@dihedral:" + dihedralName + "  @atom:" + caseJKeysI[i] + " @atom:" + caseJKeysJ[j] + " @atom:" + caseJKeysK[k] + " @atom:*\n")
            barrierEnergy = str(round(CASE_J_ENERGY / (caseJDictionaryJ[caseJKeysJ[j]] * caseJDictionaryK[caseJKeysK[k]]), 4))
            file_coeff.write("dihedral_coeff @dihedral:" + dihedralName + " " + barrierEnergy + " " + str(CASE_J_MULTIPLICITY) + " " + str(CASE_J_PHASE_SHIFT) + " 0.000\n")
            typeNames.append(dihedralName)

file_type.close()
file_coeff.close()

# Find repeated type names
len(typeNames)
repeatNames = set([x for x in typeNames if typeNames.count(x) > 1])
# Raises error if names are repeated
if len(repeatNames) != 0:
    raise Exception("Repeated type names exist in this file.")
repeatNames
