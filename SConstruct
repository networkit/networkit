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

env = Environment()

# libraries
env.Append(LIBPATH = ["", "~/workspace/gtest/lib"])
env.Append(LIBS = ["STINGER", "gtest", "log4cxx"])

env.Append(CCFLAGS = "-O0 -g3 -Wall -c -fmessage-length=0 -fopenmp -std=c++11")



# TODO: for gcc-4.6 env.Append(CCFLAGS = "-O0 -g3 -Wall -c -fmessage-length=0 -fopenmp -std=c++11")


# env.Program("EnsembleClustering-DPar", source)