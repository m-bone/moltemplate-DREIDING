# Import packages
import os
from Dreiding_label_dictionary import labelDict

os.chdir("/home/matt/Documents/Dreiding_forcefield")

# keys = list(labelDict.keys())
#
#
# acceptorKeys = [ha for ha in keys if "ha" in ha]
# donorKeys = [hd for hd in keys if "hd" in hd]

R = "2.75" # Angstroms
D = "4.0" # kcal/mol


file = open("hb_nonbond.txt", "w")

file.write("pair_coeff  @atom:*_hd  @atom:*_hd hbond/dreiding/lj @atom:H_HB i " + D + " " + R + " 4" "\n")
file.write("pair_coeff  @atom:*_hd  @atom:*_ha hbond/dreiding/lj @atom:H_HB i " + D + " " + R + " 4" "\n")

# # Donor and acceptor interactions - label = hd
# for i in range(len(donorKeys)):
#     for j in range(len(donorKeys)):
#         file.write("pair_coeff  @atom:" + donorKeys[i] + "  @atom:" + donorKeys[j] + " hbond/dreiding/lj " + "@atom:H_HB i " + D + " " + R + " 4" "\n")
#
# # Donor interactions with acceptors - label = ha
# for i in range(len(donorKeys)):
#     for j in range(len(acceptorKeys)):
#         file.write("pair_coeff  @atom:" + donorKeys[i] + "  @atom:" + acceptorKeys[j] + " hbond/dreiding/lj " + "@atom:H_HB i " + D + " " + R + " 4" "\n")

file.close()
