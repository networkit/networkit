import sys
import os
import fnmatch
from datetime import datetime

# run command on all graphs
args = sys.argv[1:]

commandTemplate = args[0]
dir = args[1]
if len(args) > 2:
    runs = int(args[2])
else:
    runs = 1
    
print("performing %d runs each" % runs)

graphFiles = []
os.chdir("/Users/forigem/Downloads/binary_networks")
for (dirpath, dirnames, filenames) in os.walk(dir):
    for name in fnmatch.filter(filenames, "*.graph"):
        path = os.path.join(dirpath, name)
        graphFiles.append(path)
        
graphFiles.sort(key=lambda s: s.lower()) # sort case-insensitively

    
commands = []
for graphFile in graphFiles:
    # outFile = "output/%s-%s.txt" % (graphFile.split("/")[-1].split(".")[0], str(datetime.now()).replace(" ", "-"))
    command = commandTemplate % {"graphFile" : graphFile}
    # command = "%s &> '%s'" % (command, outFile)
    commands.append(command)
    
print("Going to call the following %d commands:" % len(commands))
for command in commands:
    print("\t %s" % command)

called = 0
for command in commands:
    for r in range(runs):
        print("[BEGIN] %s" % command)
        os.system(command)
        called += 1
    
print("[DONE] called %d commands" % called)
