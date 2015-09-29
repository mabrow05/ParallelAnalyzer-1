#!/usr/bin/python

import os

runs = [17238,17370,17521,17925,18361,18621,18749,19232,19359,19511,19857,19899]
#runs = [19232,19359,19511,19857,19899]

sources = ["Ce139", "Sn113", "Bi207"]

runPeriod=1

if 0:
    for run in runs:
        srcPeaksByPMT = [] #holds the source peaks by PMT to be written to file
        for n in range(0,8,1):
            srcPeaksByPMT.append([])

        outfile1 = open("fits/weightedSimPeaks_runPeriod_%i.dat"%runPeriod,'w')
        outfile2 = open("fits/weightedSimPeaks_PMTbyPMT_runPeriod_%i.dat"%runPeriod,'w')
        for source in sources:
            os.system("root -l -b -q 'weightPeaks.C(%i,\"%s\")'"%(run,source))

        # to combine the peaks from a run into a source run period file
            if 0:           
                filename = "fits/%i_%s_weightedSimPeaks.dat"%(run,source)
                peakFile = open(filename,'r')
                lines1 = []
                for line in peakFile:
                    lines1.append(line)
                #print line
                peakFile.close()
                for i in range(len(lines1),0,-1):
                    outfile1.write(lines1[i-1]) # to print lower Bi first
                #print lines[i-1][:-1]
                    outfile1.write("\n")

                filename = "fits/%i_%s_weightedSimPeaks_PMTbyPMT.dat"%(run,source)
                peakFile = open(filename,'r')
                pmt=0       
                for line in peakFile:
                    peaks = line.split()
                    for i in range(len(peaks),0,-1):
                        srcPeaksByPMT[pmt].append(peaks[i-1])
                    pmt+=1
                    peakFile.close()

        for p in range(0,len(srcPeaksByPMT),1):
            for pp in range(0,len(srcPeaksByPMT[p]),1):
                outfile2.write(srcPeaksByPMT[p][pp])
                if pp<3:
                    outfile2.write(" ")
                else:
                    outfile2.write("\n")


        runPeriod+=1
        outfile1.close()
        outfile2.close()



if 1:
    srcs = ["In114E","In114W"]
    runs = [17370,17521,17925,18361,18621,18749,19232,19359,19511,19857,19899]
    for src in srcs:
        runPeriod=1
        for run in runs:
            os.system("root -l -b -q 'weightPeaks.C(%i,\"%s\")'"%(run,src))
            os.system("cp fits/%i_%s_weightedSimPeaks_PMTbyPMT.dat fits/weightedSimPeaks_PMTbyPMT_%s_runPeriod_%i.dat"%(run,src,src,runPeriod))    
            runPeriod+=1
                  
