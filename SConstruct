import os
import fnmatch

# select source
source = []

# walk source directory and find .cpp and .h
for (dirpath, dirnames, filenames) in os.walk("src"):
    for name in fnmatch.filter(filenames, "*.h"):
        source.append(os.path.join(dirpath, name))
    for name in fnmatch.filter(filenames, "*.cpp"):
        source.append(os.path.join(dirpath, name))


# exclude files matching following patterns
xpatterns = ["*ParametrizedGTest.*"]
excluded = []

for pattern in xpatterns:
	for name in fnmatch.filter(source, pattern):
		excluded.append(name)

source = [name for name in source if name not in excluded]


# set up environment (compiler flags etc.)

# 1. Debug Parallel
env = Environment()
env.Append(CCFLAGS = "-O0 -g3 -Wall -c -fmessage-length=0 -fopenmp -std=c++11")


    

# env.Program("EnsembleClustering-DPar", source)