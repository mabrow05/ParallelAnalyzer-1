import os

for run in range(16000, 20000,1):
    filename = os.getenv("REPLAY_PASS4")+"replay_pass4_%i.root"%run
    if os.path.isfile(filename):
        os.system("./rootfile_translator %i"%run)
