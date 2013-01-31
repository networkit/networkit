import os
import fnmatch

# SOURCE
source = []

# walk source directory and find ONLY .cpp files
for (dirpath, dirnames, filenames) in os.walk("src"):
    for name in fnmatch.filter(filenames, "*.cpp"):
        source.append(os.path.join(dirpath, name))


# exclude files matching following patterns
xpatterns = []
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


# TODO: extract environment-independent flags


## environment: compute

compute = Environment()
### include
compute.Append(CPPPATH = ["/home/staudt/workspace/gtest/include", \
                          "/home/staudt/workspace/STINGER/include"])
compute.Append(CCPATH = ["/home/staudt/workspace/gtest/include", \
                          "/home/staudt/workspace/STINGER/include"])
print("compute CPPPATH: %s" % compute["CPPPATH"])

### link
compute.Append(LIBS = ["STINGER", "gtest", "log4cxx"])
compute.Append(LIBPATH = ["/home/staudt/workspace/STINGER",\
                           "/home/staudt/workspace/gtest/lib"])
compute.Append(LINKFLAGS = ["-fopenmp", "-std=c++11"])

### compiler & flags
compute["CC"] = "gcc-4.7"
compute["CXX"] = "g++-4.7"




## select environment
# custom command line options
AddOption("--machine",
          dest="machine",
          type="string",
          nargs=1,
          action="store",
          help="specify the machine (environment) on which to build")


environments = {"macbook" : macbook, "compute" : compute}

try:
    env = environments[GetOption("machine")]
except:
    print("ERROR: In order to build call scons with --machine=<MACHINE> where <MACHINE> is one of: %s" % environments.keys())
    exit()


## CONFIGURATIONS

commonCFlags = ["-c", "-fmessage-length=0", "-std=c99"]
commonCppFlags = ["-std=c++11", "-Wall", "-c", "-fmessage-length=0", "-fopenmp"]

debugCppFlags = ["-O0", "-g3", "-pg"]
debugCFlags = ["-O0", "-g3"]

optimizedCppFlags = ["-O3", "-DNDEBUG"]
optimizedCFlags = ["-O3"]

# select configuration
# custom command line options
AddOption("--buildconf",
          dest="buildconf",
          type="string",
          nargs=1,
          action="store",
          help="specify the buildconfuration to build (Debug, Release)")


try:
    buildconf = GetOption("buildconf")
except:
    print("ERROR: Missing option --buildconf=<CONF>")
    exit()

# append flags

#commmon flags
env.Append(CFLAGS = commonCFlags)
env.Append(CPPFLAGS = commonCppFlags)

# buildconf flags
if buildconf == "debug":
    env.Append(CFLAGS = debugCFlags)
    env.Append(CPPFLAGS = debugCppFlags)
elif buildconf == "optimized":
    env.Append(CFLAGS = optimizedCFlags)
    env.Append(CPPFLAGS = optimizedCppFlags)
else:
    print("ERROR: invalid buildconf: %s" % buildconf)

# TARGET
env.Program("EnsembleClustering-scons-%s" % buildconf, source)