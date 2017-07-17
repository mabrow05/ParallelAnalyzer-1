#!/usr/bin/python
import sys
import csv
import os
from math import *

totfiles = 2000
numruns = 50
numfiles = totfiles/numruns

year = "2012-2013"
eastThresh = 1.19
westThresh = 1.09

jobfile = None
filestart=0
filestop=numfiles-1

with open("%s_jobs.txt"%year,"w") as jobfile:
    for i in range(0,int(numruns)):
        jobfile.write("./EfficiencyAsymm.exe %f %f %s %i %i > %s_log.txt\n"
                      %(eastThresh,westThresh,year,filestart,filestop,year))
        filestart+=numfiles
        filestop+=numfiles

os.system("parallel -P 2 < %s_jobs.txt"%(year))

