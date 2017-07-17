#!/usr/bin/python
import sys
import csv
import os
from math import *

totfiles = 10000
numruns = 100
numfiles = totfiles/numruns

field = "dip"
year = "2012-2013"

jobfile = None
filestart=0
filestop=numfiles-1

with open("%s_%s_jobs.txt"%(year,field),"w") as jobfile:
    for i in range(0,int(numruns)):
        jobfile.write("./FieldDipSystematic.exe %s %s %i %i > %s_%s_log.txt\n"
                      %(field,year,filestart,filestop,year,field))
        filestart+=numfiles
        filestop+=numfiles

os.system("parallel -j16 < %s_%s_jobs.txt"%(year,field))

