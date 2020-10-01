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
# This file generates moltemplate coefficient data for the improper angle section
# of the force field. DREIDING has a very simple implementation of improper angles
# so this script is short and doesn't use the Dreiding_label_dictionary. Improper
# angles are only specified for X_2, X_R and the C_31 united atom label.
##############################################################################

# Import packages
import os

os.chdir("/home/matt/Documents/Dreiding_forcefield")

# From the umbrella.cpp file it appears that the first atom is the middle one
#           2
#           *
#           |
#         _.*._
#       *'  0  `*
#      1         3

file = open("improper.txt", "w")

file.write("@improper:X_2  @atom:*_2* @atom:* @atom:* @atom:*\n")
file.write("@improper:X_R  @atom:*_R* @atom:* @atom:* @atom:*\n")
file.write("@improper:C_31 @atom:C_31* @atom:* @atom:* @atom:*\n")
file.write("improper_coeff @improper:X_2 40 0\n")
file.write("improper_coeff @improper:X_R 40 0\n")
file.write("improper_coeff @improper:C_31 40 54.74\n")

file.close()
