#!/usr/bin/python
import sys
import csv
import os
from math import *

totfiles = 400
numruns = 25
numfiles = totfiles/numruns

field = "flat"
year = "2011-2012"

jobfile = None
filestart=0
filestop=numfiles-1

with open("%s_%s_jobs.txt"%(year,field),"w") as jobfile:
    for i in range(0,int(numruns)):
        jobfile.write("./FieldDipSystematic_Erecon.exe %s %s %i %i > %s_%s_log.txt\n"
                      %(field,year,filestart,filestop,year,field))
        filestart+=numfiles
        filestop+=numfiles

os.system("parallel -j5 < %s_%s_jobs.txt"%(year,field))

