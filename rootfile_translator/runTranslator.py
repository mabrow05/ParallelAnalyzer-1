import os

if 0:
    for run in range(17125, 19966,1):
        filename = os.getenv("REPLAY_PASS4")+"replay_pass4_%i.root"%run
        if os.path.isfile(filename):
            os.system("./rootfile_translator %i"%run)

if 1:
    for run in range(17125, 19966,1):
        filename = os.getenv("REPLAY_PASS4")+"replay_pass4_%i.root"%run
        if os.path.isfile(filename):
            os.system("./histoAdder %i"%run)
