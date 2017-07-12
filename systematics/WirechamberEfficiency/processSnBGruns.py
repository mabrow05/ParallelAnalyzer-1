import os
import sys

bgruns = [19901,19903,19926,19928,19930,19953,19955,19957,19960]

for rn in bgruns:
    #os.system("cd ../../pedestals/; ./pedestals.exe %i"%rn)
    #os.system("cd ../../replay_pass1/; ./replay_pass1.exe %i"%rn)
    #os.system("cd ../../gain_bismuth/; ./gain_bismuth.exe %i"%rn)
    os.system("cd ../../replay_pass2/; ./replay_pass2.exe %i false"%rn)
    os.system("cd ../../replay_pass3/; ./replay_pass3.exe %i 0"%rn)
    
