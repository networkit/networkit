import os
import fnmatch

home_path = os.environ['HOME']



# SOURCE
srcDir = "src"
source = []

# walk source directory and find ONLY .cpp files
for (dirpath, dirnames, filenames) in os.walk(srcDir):
    for name in fnmatch.filter(filenames, "*.cpp"):
		source.append(os.path.join(dirpath, name))


# exclude files matching following patterns
xpatterns = ["*-X.cpp"]	# exclude executables, include later depending on target
excluded = []

for pattern in xpatterns:
	for name in fnmatch.filter(source, pattern):
		excluded.append(name)

print("excluded source files: {0}".format(excluded))
source = [name for name in source if name not in excluded]





# ENVIRONMENT

## environment: macbook
macbook = Environment()
### include
macbook.Append(CPPPATH = ["/usr/local/Cellar/gcc/4.7.2/gcc/include/c++/4.7.2", \
                          "/Users/cls/workspace/gtest/include", \
                          "/usr/local/Cellar/log4cxx/0.10.0/include"])
macbook.Append(CCPATH = ["/usr/local/Cellar/gcc/4.7.2/gcc/include/c++/4.7.2", \
                          "/Users/cls/workspace/gtest/include", \
                          "/usr/local/Cellar/log4cxx/0.10.0/include"])


### link
macbook.Append(LIBS = ["gtest", "log4cxx"])
macbook.Append(LIBPATH = ["/usr/lib/gtest", \
                            "/usr/local/Cellar/log4cxx/0.10.0/lib"])
macbook.Append(LINKFLAGS = ["-std=c++11"])

### compiler & flags
macbook["CC"] = "gcc-4.7"
macbook["CXX"] = "g++-4.7"

## environment: lappy
lappy = Environment()
### include
#myenv.Append(CPPPATH = [])
#lappy.Append(CPPPATH = ["/usr/local/Cellar/gcc/4.7.2/gcc/include/c++/4.7.2", \
#                          "/Users/cls/workspace/gtest/include", \
#                          "/usr/local/Cellar/log4cxx/0.10.0/include"])
#lappy.Append(CCPATH = ["/usr/local/Cellar/gcc/4.7.2/gcc/include/c++/4.7.2", \
#                          "/Users/cls/workspace/gtest/include", \
#                          "/usr/local/Cellar/log4cxx/0.10.0/include"])


### link
lappy.Append(LIBS = ["gtest", "log4cxx"])
lappy.Append(LIBPATH = ["/usr/src/gtest"])						# TODO: add library paths here # should work like this
lappy.Append(LINKFLAGS = ["-std=c++11"])

### compiler & flags
lappy["CC"] = "gcc-4.8"
lappy["CXX"] = "g++-4.8"


## environment: compute

compute = Environment()
### include
compute.Append(CPPPATH = [os.path.join(home_path, "workspace/gtest/include")])
compute.Append(CCPATH = [os.path.join(home_path, "workspace/gtest/include")])
print("compute CPPPATH: %s" % compute["CPPPATH"])

### link
compute.Append(LIBS = ["gtest", "log4cxx"])
compute.Append(LIBPATH = [os.path.join(home_path, "workspace/gtest")])
compute.Append(LINKFLAGS = ["-std=c++11"])

### compiler & flags
compute["CC"] = "gcc-4.7"
compute["CXX"] = "g++-4.7"

# preprocessor defines
compute.Append(CPPDEFINES=["_GNU_SOURCE"]) 


## environment: ic2.scc.kit.edu

ic2 = Environment()
### include
ic2.Append(CPPPATH = [os.path.join(home_path, "workspace/gtest/include")])
ic2.Append(CCPATH = [os.path.join(home_path, "workspace/gtest/include")])
print("compute CPPPATH: %s" % compute["CPPPATH"])

### link
ic2.Append(LIBS = ["gtest"])
ic2.Append(LIBPATH = [os.path.join(home_path, "workspace/gtest")])
ic2.Append(LINKFLAGS = ["-std=c++11"])

ic2.Append(CPPDEFINES=['NOLOG4CXX'])    # log4cxx is not available

### compiler & flags
ic2["CC"] = "/opt/gcc_4.7/bin/gcc"
ic2["CXX"] = "/opt/gcc_4.7/bin/g++"





## environment: comp_hm

comp_hm = Environment()
### include
comp_hm.Append(CPPPATH = ["/home/henningm/workspace/gtest/include", \
                          "/home/henningm/workspace/STINGER/include"])
comp_hm.Append(CCPATH = ["/home/henningm/workspace/gtest/include", \
                          "/home/henningm/workspace/STINGER/include"])
print("comp_hm CPPPATH: %s" % comp_hm["CPPPATH"])

### link
comp_hm.Append(LIBS = ["STINGER", "gtest", "log4cxx"])
comp_hm.Append(LIBPATH = ["/home/henningm/workspace/STINGER",\
                           "/home/henningm/workspace/gtest/"])
comp_hm.Append(LINKFLAGS = ["-std=c++11"])

### compiler & flags
comp_hm["CC"] = "gcc-4.7"
comp_hm["CXX"] = "g++-4.7"





## environment: mac_hm

mac_hm = Environment()
### include
mac_hm.Append(CPPPATH = ["/opt/local/include/gcc/c++/4.7.2", \
                          "/Users/Henning/Documents/workspace/gtest/include", \
                          "/opt/local/include/log4cxx/", \
                          "/Users/Henning/Documents/workspace/STINGER/include"])
mac_hm.Append(CCPATH = ["/opt/local/include/gcc/c++/4.7.2", \
                          "/Users/Henning/Documents/workspace/gtest/include", \
                          "/opt/local/include/log4cxx/", \
                          "/Users/Henning/Documents/workspace/STINGER/include"])

### link
mac_hm.Append(LIBS = ["STINGER", "gtest", "log4cxx"])
mac_hm.Append(LIBPATH = ["/Users/Henning/Documents/workspace/STINGER/OpenMP Debug",\
                           "/Users/Henning/Documents/workspace/gtest/", \
                            "/opt/local/lib/"])
mac_hm.Append(LINKFLAGS = ["-std=c++11"])
# mac_hm.Program('build/CommunityDetection.cpp')

### compiler & flags
mac_hm["CC"] = "gcc-mp-4.7"
mac_hm["CXX"] = "g++-mp-4.7"




## select environment
# custom command line options
AddOption("--machine",
          dest="machine",
          type="string",
          nargs=1,
          action="store",
          help="specify the machine (environment) on which to build")


environments = {"macbook" : macbook, "compute" : compute, "mac_hm" : mac_hm, "comp_hm" : comp_hm, "ic2": ic2, "lappy" : lappy}

try:
    env = environments[GetOption("machine")]
except:
    print("ERROR: In order to build call scons with --machine=<MACHINE> where <MACHINE> is one of: %s" % environments.keys())
    exit()


## CONFIGURATIONS

commonCFlags = ["-c", "-fmessage-length=0", "-std=c99", "-fPIC"]
commonCppFlags = ["-std=c++11", "-Wall", "-c", "-fmessage-length=0", "-fPIC"]

debugCppFlags = ["-O0", "-g3"]
debugCFlags = ["-O0", "-g3"]

optimizedCppFlags = ["-O3", "-DNDEBUG"]
optimizedCFlags = ["-O3"]

profileCppFlags = ["-O2", "-DNDEBUG", "-g", "-pg"]
profileCFlags = ["-O2", "-DNDEBUG", "-g", "-pg"]


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


# create build directory for build configuration

buildDir = ".build{0}".format(buildconf)
VariantDir(buildDir, srcDir, duplicate=0)


# modify source paths for build directory
source = [name.replace(srcDir + "/", buildDir + "/") for name in source]
print(source)


# append flags

#commmon flags
env.Append(CFLAGS = commonCFlags)
env.Append(CPPFLAGS = commonCppFlags)

# openmp yes or no
AddOption("--openmp",
          dest="openmp",
          type="string",
          nargs=1,
          action="store",
          help="-fopenmp: yes or no")

openmp = GetOption("openmp")

if (openmp == "yes") or (openmp == None): # with OpenMP by default
    env.Append(CPPFLAGS = ["-fopenmp"])
    env.Append(LINKFLAGS = ["-fopenmp"])
elif (openmp == "no"):
    env.Append(LIBS = ["pthread"])
else:
    print("ERROR: unrecognized option --openmp=%s" % openmp)
    exit()

# buildconf flags
if buildconf == "D":
    env.Append(CFLAGS = debugCFlags)
    env.Append(CPPFLAGS = debugCppFlags)
elif buildconf == "O":
    env.Append(CFLAGS = optimizedCFlags)
    env.Append(CPPFLAGS = optimizedCppFlags)
elif buildconf == "P":
	 env.Append(CFLAGS = profileCFlags)
	 env.Append(CPPFLAGS = profileCppFlags)
else:
    print("ERROR: invalid buildconf: %s" % buildconf)
    exit()


# TARGET
# openmp yes or no
AddOption("--target",
          dest="target",
          type="string",
          nargs=1,
          action="store",
          help="select target to build")


target = GetOption("target")
targetName = "NetworKit-{0}-{1}".format(target, buildconf)

if target == "CommunityDetection":
	source.append(os.path.join(srcDir, "CommunityDetection-X.cpp"))
	env.Program(targetName, source)
elif target == "DynCD":
	source.append(os.path.join(srcDir, "DynamicCommunityDetection-X.cpp"))
	env.Program(targetName, source)
elif target == "SelCD":
	source.append(os.path.join(srcDir, "SelectiveCommunityDetection-X.cpp"))
	env.Program(targetName, source)
elif target == "Core":
	# do not append executable
	env.Append(CPPDEFINES=["NOLOGGING", "NOGTEST"])
	env.Library("NetworKit-Core-{0}".format(buildconf), source)
else:
	print("ERROR: unknown target: %" % target)
	exit()





# TODO: make unit tests a separate target
