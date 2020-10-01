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
