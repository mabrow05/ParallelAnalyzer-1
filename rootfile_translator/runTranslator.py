import os

for run in range(17125, 19966,1):
    filename = os.getenv("REPLAY_PASS4")+"replay_pass4_%i.root"%run
    if os.path.isfile(filename):
        os.system("./rootfile_translator %i"%run)
