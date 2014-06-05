import os
import fnmatch
import ConfigParser
import distutils.sysconfig

home_path = os.environ['HOME']

# SOURCE files (including executable) will be gathered here
srcDir = "src/cpp"
def getSourceFiles(target, optimize):
	source = []

	# walk source directory and find ONLY .cpp files
	for (dirpath, dirnames, filenames) in os.walk(srcDir):
		for name in fnmatch.filter(filenames, "*.cpp"):
			source.append(os.path.join(dirpath, name))

	# exclude files depending on target, executables will be addes later
	if (target not in ["Tests"]):
		# exclude files matching following patterns
		xpatterns = ["*-X.cpp","*Unittests.cpp","*GTest.cpp","*Benchmark.cpp"]
		excluded = []
	else:
		# exclude files matching following patterns
		# only the target "Test" requires Benchmark and GTest files
		xpatterns = ["*-X.cpp","*Unittests.cpp"]
		excluded = []

	for pattern in xpatterns:
		for name in fnmatch.filter(source, pattern):
			excluded.append(name)

	#print("excluded source files: {0}".format(excluded))
	source = [name for name in source if name not in excluded]

	# add executable
	if target == "Tests":
		source.append(os.path.join(srcDir, "Unittests-X.cpp"))
	elif target == "Core":
		pass # no executable
	# elif target == "Python":
	# 	# add SWIG interface files
	# 	for (dirpath, dirnames, filenames) in os.walk("src/swig"):
	# 		for name in fnmatch.filter(filenames, "*.i"):
	# 			source.append(os.path.join(dirpath, name))
	# else:
	# 	print("Unknown target: {0}".format(target))
	# 	Exit(1)

	# create build directory for build configuration
	buildDir = ".build{0}".format(optimize)
	VariantDir(buildDir, srcDir, duplicate=0)

	# modify source paths for build directory
	source = [name.replace(srcDir + "/", buildDir + "/") for name in source]
	#print(source)
	return source


# ENVIRONMENT

## read environment settings from configuration file

env = Environment()
confPath = "build.conf"
if not os.path.isfile(confPath):
	print("The configuration file `build.conf` does not exist. You need to create it.")
	print("Use the file build.conf.example to create your build.conf")
	Exit(1)

conf = ConfigParser.ConfigParser()
conf.read([confPath])     # read the configuration file

## compiler
cppComp = conf.get("compiler", "cpp", "gcc")
defines = conf.get("compiler", "defines", [])		# defines are optional
if defines is not []:
    defines = defines.split(",")

## includes
stdInclude = conf.get("includes", "std", "")      # includes for the standard library - may not be needed
gtestInclude = conf.get("includes", "gtest")
if conf.has_option("includes", "tbb"):
	tbbInclude = conf.get("includes", "tbb", "")
else:
	tbbInclude = ""

## libraries
gtestLib = conf.get("libraries", "gtest")
if conf.has_option("libraries", "tbb"):
	tbbLib = conf.get("libraries", "tbb", "")
else:
	tbbLib = ""

pythonInclude = "-I/usr/local/Cellar/python3/3.4.1/Frameworks/Python.framework/Versions/3.4/include/python3.4m"
pythonInclude = "-I/usr/include/python2.7"

env["CC"] = cppComp
env["CXX"] = cppComp

env.Append(CPPDEFINES = defines)
env.Append(CPPPATH = [stdInclude, gtestInclude, tbbInclude, pythonInclude])
env.Append(LIBS = ["gtest"])
env.Append(LIBPATH = [gtestLib, tbbLib])
env.Append(LINKFLAGS = ["-std=c++11"])

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
AddOption("--optimize",
          dest="optimize",
          type="string",
          nargs=1,
          action="store",
          help="specify the optimization level to build: D(ebug), O(ptimize), P(rofile)")

AddOption("--sanitize",
          dest="sanitize",
          type="string",
          nargs=1,
          action="store",
          help="switch on address sanitizer")


try:
    optimize = GetOption("optimize")
except:
    print("ERROR: Missing option --optimize=<LEVEL>")
    exit()

sanitize = None
try:
	sanitize = GetOption("sanitize")
except:
	pass
	


# create build directory for build configuration
# modify source paths for build directory
# moved to getSourceFiles()

# append flags

#commmon flags
env.Append(CFLAGS = commonCFlags)
env.Append(CPPFLAGS = commonCppFlags)

# logging yes or no
AddOption("--logging",
          dest="logging",
          type="string",
          nargs=1,
          action="store",
          help="enable logging: yes or no")

logging = GetOption("logging")

if logging == "no":
	env.Append(CPPDEFINES=["NOLOGGING"]) # logging is enabled by default
	print("INFO: Logging is now disabled")
elif (logging != "yes") and (logging != None):
	print("INFO: unrecognized option --logging=%s" % logging)
	print("Logging is enabled by default")

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

# optimize flags
if optimize == "D":
    env.Append(CFLAGS = debugCFlags)
    env.Append(CPPFLAGS = debugCppFlags)
elif optimize == "O":
    env.Append(CFLAGS = optimizedCFlags)
    env.Append(CPPFLAGS = optimizedCppFlags)
elif optimize == "P":
	 env.Append(CFLAGS = profileCFlags)
	 env.Append(CPPFLAGS = profileCppFlags)
else:
    print("ERROR: invalid optimize: %s" % optimize)
    exit()

# sanitize
if sanitize:
	if sanitize == "address":
		env.Append(CPPFLAGS = ["-fsanitize=address"])
		env.Append(LINKFLAGS = ["-fsanitize=address"])
	else:
		print("ERROR: invalid sanitize option")
		exit()


# TARGET
AddOption("--target",
          dest="target",
          type="string",
          nargs=1,
          action="store",
          help="select target to build")


target = GetOption("target")
availableTargets = ["Core", "Tests", "Python"]
if target in availableTargets:
	source = getSourceFiles(target,optimize)
	targetName = "NetworKit-{0}-{1}".format(target, optimize)
	if target == "Core":
		env.Library(targetName, source)
	elif target == "Tests":
		env.Program(targetName, source)
	elif target == "Python":
		import distutils.sysconfig
		for (dirpath, dirnames, filenames) in os.walk("src/swig"):
			for name in fnmatch.filter(filenames, "*.i"):
				source.append(os.path.join(dirpath, name))
		env.Append(SWIGFLAGS=['-c++', '-python'])
		env.Append(SHLIBPREFIX="")
		env.Append(LINKFLAGS = ["-lpython"])
		env.SharedLibrary('_NetworKit.so', source)
		Command(target = "src/swig/_NetworKit.so",
	        source = "lib_NetworKit.so",
	        action = [Delete("$TARGET"), Move("$TARGET", "$SOURCE")])

   	elif target == "Python2":
		env.Append(SWIGFLAGS = ["-c++", "-python"])
		env.Append(tools = ['default', 'swig'])
		env.SharedLibrary("src/swig/_NetworKit", source)

		for (dirpath, dirnames, filenames) in os.walk("src/swig"):
			for name in fnmatch.filter(filenames, "*.i"):
				source.append(os.path.join(dirpath, name))
		print(source)

		import distutils.sysconfig
		# env.Append(CPPPATH=[distutils.sysconfig.get_python_inc()])
		# env.Append(SHLIBPREFIX="")
		env.SharedLibrary('_example.so', source)

		# Command(target = "src/swig/_NetworKit.so",
	        # source = "src/swig/lib_NetworKit.so",
	        # action = [Delete("$TARGET"), Move("$TARGET", "$SOURCE")])	
		# Command(target = "src/swig/_NetworKit.dylib",
	 #        source = "src/swig/lib_NetworKit.dylib",
	 #        action = [Delete("$TARGET"), Move("$TARGET", "$SOURCE")])
else:
	print("ERROR: unknown target: {0}".format(target))
	exit()
