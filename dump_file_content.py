#script to dump the content of a file

import uproot
import sys


file_path = "/pnfs/psi.ch/cms/trivcat/store/user/cazzanig/PFNano_s-channel_mMed-1500_mDark-20_rinv-0.3_alpha-peak_13TeV-pythia8_n-1000_part-1.root"

file = uproot.open(file_path)

print(file.keys())

print(file["mmtree/Events"].show())

#dump the content in a .txt file
with open("dump.txt", "w") as f:
    f.write(file["mmtree/Events"].show())