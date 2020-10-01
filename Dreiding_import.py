import os
import pandas as pd
import numpy as np

#os.chdir("/home/matt/Documents/Dreiding_forcefield")

input = pd.read_csv("phencharge.lt",  header = None)
input

atomDataStart = input[0].str.contains('Atoms')
atomDataStart = atomDataStart[atomDataStart].index


atomDataEnd = input[0].str.contains('}')
atomDataEnd = atomDataEnd[atomDataEnd].index
atomDataEnd = atomDataEnd[atomDataEnd > atomDataStart[0]][0]

atomData = input[atomDataStart[0] + 1: atomDataEnd]
atomData = atomData.replace("\\t", "", regex = True)
atomData.to_csv("data.txt", index = False, header = False)

reinput = pd.read_csv("data.txt", sep = " ",  header = None)
reinput = reinput.dropna(1)

data = reinput.iloc[:, 2:4]
