#!/usr/bin/python
import sys
import csv
import os
from math import *

totfiles = 2000
numruns = 50
numfiles = totfiles/numruns

year = "2011-2012"
eastThresh = 1.019#0.2#
westThresh = 1.019#0.2#

jobfile = None
filestart=0
filestop=numfiles-1

with open("%s_E%0.2f_W%0.2f_jobs.txt"%(year,eastThresh,westThresh),"w") as jobfile:
    for i in range(0,int(numruns)):
        jobfile.write("./EfficiencyAsymm.exe %f %f %s %i %i > %s_log.txt\n"
                      %(eastThresh,westThresh,year,filestart,filestop,year))
        filestart+=numfiles
        filestop+=numfiles

os.system("parallel -P 6 < %s_E%0.2f_W%0.2f_jobs.txt"%(year,eastThresh,westThresh))

