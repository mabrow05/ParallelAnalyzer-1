import os
import sys
from math import *

#load uncorrected asymm by bin
enlow = []
uncorr = []
with open("BScorrs/UNBLINDED_DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_0-59_BinByBin.txt","r") as ifile:
    for line in ifile:
        l=line.split("\t")
        uncorr.append(float(l[1]))
        enlow.append(float(l[0]))

#load deltaBS0
deltaBS0 = []
with open("BScorrs/UNBLINDED_DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_0-59_DeltaBS_0_BinByBin.txt","r") as ifile:
    for line in ifile:
        l=line.split("\t")
        deltaBS0.append(float(l[1]))


#load deltaBS1
deltaBS1 = []
with open("BScorrs/UNBLINDED_DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_0-59_DeltaBS_1_BinByBin.txt","r") as ifile:
    for line in ifile:
        l=line.split("\t")
        deltaBS1.append(float(l[1]))


#load deltaBS2
deltaBS2 = []
with open("BScorrs/UNBLINDED_DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_0-59_DeltaBS_2_BinByBin.txt","r") as ifile:
    for line in ifile:
        l=line.split("\t")
        deltaBS2.append(float(l[1]))


#load deltaBS3
deltaBS3 = []
with open("BScorrs/UNBLINDED_DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_0-59_DeltaBS_3_BinByBin.txt","r") as ifile:
    for line in ifile:
        l=line.split("\t")
        deltaBS3.append(float(l[1]))


#write corr files
with open("deltaBS0_2011-2012.txt","w") as f0,open("deltaBS1_2011-2012.txt","w") as f1,open("deltaBS2_2011-2012.txt","w") as f2,open("deltaBS3_2011-2012.txt","w") as f3:
    for i in range(0,len(uncorr)):
        en = enlow[i]
        f0.write("%f\t%0.10f\t%0.10f\n"%(en,deltaBS0[i]/uncorr[i]-1. if uncorr[i]!=0. else 0,fabs(deltaBS0[i]/uncorr[i]-1.)*0.25 if uncorr[i]!=0. else 0))
        f1.write("%f\t%0.10f\t%0.10f\n"%(en,deltaBS1[i]/deltaBS0[i]-1. if deltaBS0[i]!=0. else 0,fabs(deltaBS1[i]/deltaBS0[i]-1.)*0.25 if deltaBS0[i]!=0. else 0))
        f2.write("%f\t%0.10f\t%0.10f\n"%(en,deltaBS2[i]/deltaBS1[i]-1. if deltaBS1[i]!=0. else 0,fabs(deltaBS2[i]/deltaBS1[i]-1.)*0.25 if deltaBS1[i]!=0. else 0))
        f3.write("%f\t%0.10f\t%0.10f\n"%(en,deltaBS3[i]/deltaBS2[i]-1. if deltaBS2[i]!=0. else 0,fabs(deltaBS3[i]/deltaBS2[i]-1.)*0.25 if deltaBS2[i]!=0. else 0))




########################################################
# 2012-2013

#load uncorrected asymm by bin
enlow = []
uncorr = []
with open("BScorrs/UNBLINDED_DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_60-121_BinByBin.txt","r") as ifile:
    for line in ifile:
        l=line.split("\t")
        uncorr.append(float(l[1]))
        enlow.append(float(l[0]))

#load deltaBS0
deltaBS0 = []
with open("BScorrs/UNBLINDED_DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_60-121_DeltaBS_0_BinByBin.txt","r") as ifile:
    for line in ifile:
        l=line.split("\t")
        deltaBS0.append(float(l[1]))


#load deltaBS1
deltaBS1 = []
with open("BScorrs/UNBLINDED_DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_60-121_DeltaBS_1_BinByBin.txt","r") as ifile:
    for line in ifile:
        l=line.split("\t")
        deltaBS1.append(float(l[1]))


#load deltaBS2
deltaBS2 = []
with open("BScorrs/UNBLINDED_DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_60-121_DeltaBS_2_BinByBin.txt","r") as ifile:
    for line in ifile:
        l=line.split("\t")
        deltaBS2.append(float(l[1]))


#load deltaBS3
deltaBS3 = []
with open("BScorrs/UNBLINDED_DeltaTheoryOnly_OctetAsymmetries_AnaChC_Octets_60-121_DeltaBS_3_BinByBin.txt","r") as ifile:
    for line in ifile:
        l=line.split("\t")
        deltaBS3.append(float(l[1]))


#write corr files
with open("deltaBS0_2012-2013.txt","w") as f0,open("deltaBS1_2012-2013.txt","w") as f1,open("deltaBS2_2012-2013.txt","w") as f2,open("deltaBS3_2012-2013.txt","w") as f3:
    for i in range(0,len(uncorr)):
        en = enlow[i]
        f0.write("%f\t%0.10f\t%0.10f\n"%(en,deltaBS0[i]/uncorr[i]-1. if uncorr[i]!=0. else 0,fabs(deltaBS0[i]/uncorr[i]-1.)*0.25 if uncorr[i]!=0. else 0))
        f1.write("%f\t%0.10f\t%0.10f\n"%(en,deltaBS1[i]/deltaBS0[i]-1. if deltaBS0[i]!=0. else 0,fabs(deltaBS1[i]/deltaBS0[i]-1.)*0.25 if deltaBS0[i]!=0. else 0))
        f2.write("%f\t%0.10f\t%0.10f\n"%(en,deltaBS2[i]/deltaBS1[i]-1. if deltaBS1[i]!=0. else 0,fabs(deltaBS2[i]/deltaBS1[i]-1.)*0.25 if deltaBS1[i]!=0. else 0))
        f3.write("%f\t%0.10f\t%0.10f\n"%(en,deltaBS3[i]/deltaBS2[i]-1. if deltaBS2[i]!=0. else 0,fabs(deltaBS3[i]/deltaBS2[i]-1.)*0.25 if deltaBS2[i]!=0. else 0))

