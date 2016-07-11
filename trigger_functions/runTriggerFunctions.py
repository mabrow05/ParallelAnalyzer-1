#!/usr/bin/python

import os
import sys

xenonList = [4,5,7]
octetRange = [40, 59]


#for index in range(octetRange[0], octetRange[1], 1):
 #   os.system("./findADCthreshold.exe octet %i"%index)

for index in xenonList:
    os.system("./findADCthreshold.exe xenon %i"%index)
