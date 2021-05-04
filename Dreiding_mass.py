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
# This file generates moltemplate label and mass data. It simply combines the
# force field label with the atomic mass of the atom in the label. Included is a
# loop to output a list of atom labels as a text file.
##############################################################################

# Import packages
import os
from Dreiding_label_dictionary import labelDict

os.chdir("/home/matt/Documents/XP_Project/Dreiding_forcefield")

#List of keys from dictionary
keys = list(labelDict.keys())

file = open("masses.txt", "w")

#Write labels for documentation - not needed everytime, built for LaTex output
# for k in range(len(keys)):
#   file.write("\item " + keys[k].replace("_", "\_") + "\n")

repeats = 3

# Loop through all dictionary values and the number of repeats specified
for h in range(len(labelDict)):
    for i in range(1, repeats +1):
        file.write("@atom:" + str(i) + keys[h][1:] + "\t" + str(labelDict[keys[h]][0]) + "\n")

file.close()
