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
# This file generates moltemplate coefficient data for the hydrogen bonding
# component of the force field. LAMMPS uses hbond/dreiding style and requires
# specialised forcefield coefficients to faciliate it. The use of the hd/ha
# flags and wildcards (*) means this file is very short.
##############################################################################

# Import packages
import os

os.chdir("/home/matt/Documents/UsefulStore/XP_Project/Dreiding_forcefield/")

R = "2.75" # Angstroms
D = "4.0" # kcal/mol

file = open("hb_nonbond.txt", "w")

file.write("pair_coeff  @atom:*_hd  @atom:*_hd hbond/dreiding/lj @atom:H_HB i " + D + " " + R + " 4" "\n")
file.write("pair_coeff  @atom:*_hd  @atom:*_hd hbond/dreiding/lj @atom:H_HB j " + D + " " + R + " 4" "\n")
file.write("pair_coeff  @atom:*_hd  @atom:*_ha hbond/dreiding/lj @atom:H_HB i " + D + " " + R + " 4" "\n")

file.close()
