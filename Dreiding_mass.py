# Import packages
import os
from Dreiding_label_dictionary import labelDict

os.chdir("/home/matt/Documents/Dreiding_forcefield")

#List of keys from dictionary
keys = list(labelDict.keys())

# Number of labels with sequential label numbers. Result is repeats-1 labels
repeats = 1

file = open("masses.txt", "w")

#Write labels for documentation - not needed everytime, built for LaTex output
# for k in range(len(keys)):
#   file.write("\item " + keys[k].replace("_", "\_") + "\n")

# Loop through all dictionary values and the number of repeats specified
for h in range(len(labelDict)):
    for i in range(repeats):

        # Keeps most simple label free if user doesn't want label numbers
        if i == 0:
            file.write("@atom:" + keys[h] + "\t" + str(labelDict[keys[h]][0]) + "\n")
            continue

        # Inserts number after label for users to specify different atom types
        #file.write("@atom:" + keys[h] + "_" + str(i) + "\t" + str(labelDict[keys[h]])+ "\n")

file.close()
