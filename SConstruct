import os
import subprocess
import fnmatch
import ConfigParser

home_path = os.environ['HOME']

def checkStd(compiler):
	sample = open("sample.cpp", "w")
	sample.write("""
	#include <iostream>

	[[deprecated("use the function body directly instead of wrapping it in a function.")]]
	void helloWorld() {
		std::cout << "Hello world" << std::endl;
	}

	int main (int argc, char *argv[]) {
		helloWorld();
		return 0;
	}""")
	sample.close()
	FNULL = open(os.devnull, 'w')
	if subprocess.call([compiler,"-o","test_build","-std=c++14","sample.cpp"],stdout=FNULL,stderr=FNULL) == 0:
		stdflag = "c++14"
	elif subprocess.call([compiler,"-o","test_build","-std=c++11","sample.cpp"],stdout=FNULL,stderr=FNULL) == 0:
		stdflag = "c++11"
	else:
		# possibility to print warning/error
		# assume c++11
		stdflag = "c++11"
	# clean up
	FNULL.close()
	os.remove("sample.cpp")
	try:
		os.remove("test_build")
	except:
		pass
	return stdflag


# SOURCE files (including executable) will be gathered here
srcDir = "networkit/cpp"
def getSourceFiles(target, optimize):
	source = []

	# walk source directory and find ONLY .cpp files
	for (dirpath, dirnames, filenames) in os.walk(srcDir):
	    for name in fnmatch.filter(filenames, "*.cpp"):
			source.append(os.path.join(dirpath, name))

	# exclude files depending on target, executables will be addes later
	xpatterns = ["*-X.cpp"]
	excluded = []

	# only the target "Test" requires Benchmark and GTest files
	if (target not in ["Tests"]):
		# exclude files matching following patterns
		xpatterns += ["*GTest.cpp","*Benchmark.cpp"]

	for pattern in xpatterns:
		for name in fnmatch.filter(source, pattern):
			excluded.append(name)

	#print("excluded source files: {0}".format(excluded))
	source = [name for name in source if name not in excluded]

	# add executable
	if target == "Tests":
		source.append(os.path.join(srcDir, "Unittests-X.cpp"))
	elif target in ["Core","Lib"]:
		pass # no executable
	else:
		print("Unknown target: {0}".format(target))
		Exit(1)

	# create build directory for build configuration
	buildDir = ".build{0}".format(optimize)
	VariantDir(buildDir, srcDir, duplicate=0)

	# modify source paths for build directory
	source = [name.replace(srcDir + "/", buildDir + "/") for name in source]
	#print(source)
	return source


AddOption("--compiler",
	dest="compiler",
	type="string",
	nargs=1,
	action="store",
	help="used to pass gcc version from setup.py to SConstruct")

AddOption("--std",
	dest="std",
	type="string",
	nargs=1,
	action="store",
	help="used to pass std flag from setup.py to SConstruct")

# ENVIRONMENT

## read environment settings from configuration file

env = Environment()
compiler = GetOption("compiler")
stdflag = GetOption("std")

if not os.path.isfile("build.conf") and not compiler == None:
	#print("{0} has been passed via command line".format(compiler))
	env["CC"] = compiler
	env["CXX"] = compiler
else:
	confPath = "build.conf"
	if not os.path.isfile(confPath):
		print("The configuration file `build.conf` does not exist. You need to create it.")
		print("Use the file build.conf.example to create your build.conf")
		Exit(1)

	conf = ConfigParser.ConfigParser()
	conf.read([confPath])     # read the configuration file

	## compiler
	if compiler is None:
		cppComp = conf.get("compiler", "cpp", "gcc")
	else:
		cppComp = compiler
	defines = conf.get("compiler", "defines", [])		# defines are optional
	if defines is not []:
		defines = defines.split(",")


	## C++14 support
	if stdflag is None:
		try:
			stdflag = conf.get("compiler", "std14")
		except:
			pass
	if stdflag is None or len(stdflag) == 0:
		# do test compile
		stdflag = checkStd(cppComp)
		# and store it in the configuration
		conf.set("compiler","std14", stdflag)

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

	env["CC"] = cppComp
	env["CXX"] = cppComp

	env.Append(CPPDEFINES=defines)
	env.Append(CPPPATH = [stdInclude, gtestInclude, tbbInclude])
	env.Append(LIBS = ["gtest"])
	env.Append(LIBPATH = [gtestLib, tbbLib])

	with open(confPath, "w") as f:
		conf.write(f)

env.Append(LINKFLAGS = ["-std={}".format(stdflag)])

## CONFIGURATIONS

commonCFlags = ["-c", "-fmessage-length=0", "-std=c99", "-fPIC"]
commonCppFlags = ["-std={}".format(stdflag), "-Wall", "-c", "-fmessage-length=0", "-fPIC"]

debugCppFlags = ["-O0", "-g3", "-DLOG_LEVEL=LOG_LEVEL_TRACE"]
debugCFlags = ["-O0", "-g3"]

optimizedCppFlags = ["-O3", "-DNDEBUG", "-DLOG_LEVEL=LOG_LEVEL_INFO"]
optimizedCFlags = ["-O3"]

profileCppFlags = ["-O2", "-DNDEBUG", "-g", "-pg", "-DLOG_LEVEL=LOG_LEVEL_INFO"]
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
if optimize == "Dbg":
    env.Append(CFLAGS = debugCFlags)
    env.Append(CPPFLAGS = debugCppFlags)
elif optimize == "Opt":
    env.Append(CFLAGS = optimizedCFlags)
    env.Append(CPPFLAGS = optimizedCppFlags)
elif optimize == "Pro":
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
availableTargets = ["Lib","Core","Tests"]
if target in availableTargets:
	source = getSourceFiles(target,optimize)
	targetName = "NetworKit-{0}-{1}".format(target, optimize)
	if target in ["Core","Lib"]:
		# do not append executable
		# env.Append(CPPDEFINES=["NOLOGGING"])
		env.Library("NetworKit-Core-{0}".format(optimize), source)
		if target == "Lib":
			libFileToLink = "libNetworKit-Core-{0}.a".format(optimize)
			libFileTarget = "libNetworKit.a"
			if os.path.lexists(libFileTarget):
				os.remove(libFileTarget)
			os.symlink(libFileToLink,libFileTarget)
			# SCons does not support python 3 yet...
			#os.symlink("src/cpp","NetworKit",True)
			# to support case insensitive file systems
			# place the symlink for the include path in the folder include
			if os.path.isdir("include"):
				try:
					os.remove("include/NetworKit")
				except:
					pass
				os.rmdir("include")
			os.mkdir("include")
			os.chdir("include")
			subprocess.call(["ln","-s","../networkit/cpp","NetworKit"])
			os.chdir("../")

	else:
		env.Program(targetName, source)
else:
	print("ERROR: unknown target: {0}".format(target))
	exit(1)
