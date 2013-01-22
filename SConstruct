import os
import fnmatch

# SOURCE
source = []

# walk source directory and find ONLY .cpp files
for (dirpath, dirnames, filenames) in os.walk("src"):
    for name in fnmatch.filter(filenames, "*.cpp"):
        source.append(os.path.join(dirpath, name))


# exclude files matching following patterns
xpatterns = ["*ParametrizedGTest.*"]
excluded = []

for pattern in xpatterns:
	for name in fnmatch.filter(source, pattern):
		excluded.append(name)

source = [name for name in source if name not in excluded]


# ENVIRONMENT

## environment: macbook
macbook = Environment()
### include
macbook.Append(CPPPATH = ["/usr/local/Cellar/gcc/4.7.2/gcc/include/c++/4.7.2", \
                          "/Users/cls/workspace/gtest/include", \
                          "/usr/local/Cellar/log4cxx/0.10.0/include", \
                          "/Users/cls/workspace/STINGER/include"])
macbook.Append(CCPATH = ["/usr/local/Cellar/gcc/4.7.2/gcc/include/c++/4.7.2", \
                          "/Users/cls/workspace/gtest/include", \
                          "/usr/local/Cellar/log4cxx/0.10.0/include", \
                          "/Users/cls/workspace/STINGER/include"])


### link
macbook.Append(LIBS = ["STINGER", "gtest", "log4cxx"])
macbook.Append(LIBPATH = ["/Users/cls/workspace/STINGER/OpenMP Debug",\
                           "/Users/cls/workspace/gtest/lib", \
                            "/usr/local/Cellar/log4cxx/0.10.0/lib"])
macbook.Append(LINKFLAGS = ["-fopenmp", "-std=c++11"])

### compiler & flags
macbook["CC"] = "gcc-4.7"
macbook["CXX"] = "g++-4.7"
macbook.Append(CFLAGS = ["-c", "-fmessage-length=0", "-std=c99"])
macbook.Append(CPPFLAGS = ["-std=c++11", "-O0", "-g3", "-Wall", "-c", "-fmessage-length=0", "-g", "-pg", "-fopenmp"])


# TODO: extract environment-independent flags


## environment: compute11

compute11 = Environment()
### include
compute11.Append(CPPPATH = ["/home/staudt/workspace/gtest/include", \
                          "/home/staudt/workspace/STINGER/include"])
compute11.Append(CCPATH = ["/home/staudt/workspace/gtest/include", \
                          "/home/staudt/workspace/STINGER/include"])
print("compute11 CPPPATH: %s" % compute11["CPPPATH"])

### link
compute11.Append(LIBS = ["STINGER", "gtest", "log4cxx"])
compute11.Append(LIBPATH = ["/home/staudt/workspace/STINGER",\
                           "/home/staudt/workspace/gtest/lib"])
compute11.Append(LINKFLAGS = ["-fopenmp", "-std=c++11"])

### compiler & flags
compute11["CC"] = "gcc-4.7"
compute11["CXX"] = "g++-4.7"
compute11.Append(CFLAGS = ["-c", "-fmessage-length=0", "-std=c99"])
compute11.Append(CPPFLAGS = ["-std=c++11", "-O0", "-g3", "-Wall", "-c", "-fmessage-length=0", "-g", "-pg", "-fopenmp"])


# TODO: for gcc-4.6 env.Append(CCFLAGS = "-O0 -g3 -Wall -c -fmessage-length=0 -fopenmp -std=c++11")


## select environment
# custom command line options
AddOption("--machine",
          dest="machine",
          type="string",
          nargs=1,
          action="store",
          help="specify the machine (environment) on which to build")


environments = {"macbook" : macbook, "compute11" : compute11}

try:
    env = environments[GetOption("machine")]
except:
    print("ERROR: In order to build call scons with --machine=<MACHINE> where <MACHINE> is one of: %s" % environments.keys())
    exit()


# TARGET
env.Program("EnsembleClustering-DPar-scons", source)