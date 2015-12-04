#!/usr/bin/python

#produce the files which hold the number of sources and what type for each run

import os

with open("UCNA Run Log.txt") as log:
    srcs = []
    for line in log:      
        if (line[:8]=="@sources"):
            print "..."
            srcs = []
            words = line.split()# list(enumerate(line.split()))
            for r in words:
                if r[:1]!="@" and len(r)>1:# and r[len(r)-1]!="\n":
                    srcs.append(r)
                
                         
            
        if (line[7:17]=="SourcesCal" and int(line[1:7])>17000):
            if len(srcs)>0:
                run =  int(line[1:7])
                print run,srcs
                filename = "source_list_MB/source_list_%i.dat"%run
                outfile = open(filename,"w")
                outfile.write("%i\n"%len(srcs))
                for src in srcs:
                    outfile.write("%s\n"%src[:2])
                outfile.close()

