import os

runs = [17238]#,17370,17521,17925,18361,18621,18749,19232,19359,19511,19857,19899]
sources = ["Ce139", "Sn113", "Bi207"]

for run in runs:
    for source in sources:
        os.system("root -l -b -q 'weightPeaks.C(%i,\"%s\")'"%(run,source))
        
