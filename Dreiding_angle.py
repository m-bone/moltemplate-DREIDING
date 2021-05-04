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
# This file generates moltemplate coefficient data for the bond angle section of
# the force field. DREIDING has a very simple implementation of bond angles so
# this script is short and doesn't use the Dreiding_label_dictionary. The
# atomDictionary is the list of possible labels and associated angles. In
# addition, X_2, X_1 and X_R labels are added with wildcards as they all share
# the same angles (120.0, 180.0, and 120.0 respectively). The angle constant is
# fixed at 50.0 for all angle types.
##############################################################################

# Import packages
import os

os.chdir("/home/matt/Documents/XP_Project/Dreiding_forcefield")

# Dictionary of bond angles (in deg) according to DREIDING
atomDictionary = {
    "*H": 180.0,
    "*H_HB": 180.0,
    "*H_B": 90.0,
    "*B_3": 109.471,
    "*C_3": 109.471,
    "*N_3": 106.7,
    "*O_3": 104.51,
    "*F": 180.0,
    "*Al_3": 109.471,
    "*Si_3": 109.471,
    "*P_3": 93.3,
    "*S_3": 92.1,
    "*Cl": 180.0,
    "*Ga_3": 109.471,
    "*Ge_3": 109.471,
    "*As_3": 92.1,
    "*Se_3": 90.6,
    "*Br": 180.0,
    "*In_3": 109.471,
    "*Sn_3": 109.471,
    "*Sb_3": 91.6,
    "*Te_3": 90.3,
    "*I": 180.0,
    "*Na": 90.0,
    "*Ca": 90.0,
    "*Fe": 90.0,
    "*Zn": 109.471
}

#List of keys from dictionary
keys = list(atomDictionary.keys())

# Create text file
file = open("angle.txt", "w")

# Angles by type
for x in range(len(atomDictionary)):
    file.write("@angle:" + keys[x][1:] + "  @atom:* @atom:" + keys[x] + "* @atom:*\n")

# Write specific X_2, X_1 and X_R types as they all have the same angles
file.write("@angle:_2  @atom:* @atom:*_2* @atom:*\n")
file.write("@angle:_1  @atom:* @atom:*_1* @atom:*\n")
file.write("@angle:_R  @atom:* @atom:*_R* @atom:*\n")


# Angle coeffs
for y in range(len(atomDictionary)):
    file.write("angle_coeff @angle:" + keys[y][1:] + " 50.0 " + str(atomDictionary[keys[y]]) + "\n")

# Write specific X_2, X_1 and X_R types as they all have the same angles
file.write("angle_coeff @angle:_2 50.0 120.0 \n")
file.write("angle_coeff @angle:_1 50.0 180.0 \n")
file.write("angle_coeff @angle:_R 50.0 120.0 \n")

file.close()
