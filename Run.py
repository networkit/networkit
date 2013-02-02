import sys
import os
import fnmatch
from datetime import datetime

# run command on all graphs
args = sys.argv[1:]

commandTemplate = args[0]
dir = args[1]

graphFiles = []

for (dirpath, dirnames, filenames) in os.walk(dir):
    print(filenames)
    for name in fnmatch.filter(filenames, "*.graph"):
        print(name)
        path = os.path.join(dirpath, name)
        print(path)
        graphFiles.append(path)
        
print("batch: %s" % graphFiles)
        
for graphFile in graphFiles:
    outFile = "output/%s-%s.txt" % (graphFile.split("/")[-1].split(".")[0], str(datetime.now()).replace(" ", "-"))
    command = commandTemplate % {"graphFile" : graphFile}
    command = "%s > '%s'" % (command, outFile)
    os.system(command)