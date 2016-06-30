import sys

##################################
# check Python version
##################################

if sys.version_info.major < 3:
	print("ERROR: NetworKit requires Python 3.")
	sys.exit(1)

import version
from setup_util import *
if "setuptools" not in sys.modules:
	from ez_setup import use_setuptools
	# in case setuptools is not installed
	use_setuptools()

from setuptools import setup, Extension, find_packages
from setuptools.command.test import test as TestCommand
from setuptools.command.build_ext import build_ext as SetuptoolsBuildExtCmd
from setuptools.command.install import install as InstallCmd
# from setuptools.command.clean import clean as CleanCmd
import unittest

abort_installation = False
errorMessages = []
warnMessages = []
try:
	import Cython
	from Cython.Build import cythonize
	from Cython.Distutils import build_ext as CythonBuildExtCmd
	from distutils.version import LooseVersion
	if LooseVersion(Cython.__version__) >= LooseVersion('0.21'):
		cython_available = True
	else:
		cython_available = False
		#print("Cython version too old, please update")
except:
	# import so that the deriving class can still be there
	from setuptools.command.build_ext import build_ext as CythonBuildExtCmd
	abort_installation = False
	cython_available = False
	# errorMessages.append("ERROR: Cython not installed. Please install Cython and rerun")

import multiprocessing
import os
import shutil

import subprocess

from argparse import ArgumentParser

##################################
# check whether SCons is available
##################################
scons_available = None
try:
	if shutil.which("scons") is None:
		# errorMessages.append("ERROR: Build system SCons is not installed. Please install and rerun")
		# abort_installation = True
		scons_available = False
	else:
		scons_available = True
except:
	print("WARNING: unable to check whether build system SCons is installed")
	scons_available = False

#
if sys.platform == 'Windows' and not scons_available:
	abort_installation = True
	errorMessages.append("ERROR: Build system SCons is not installed. Please install and rerun")


# compiler candidates
# this list serves as a fallback when neither $CXX is set nor build.conf exists
candidates = ["g++", "g++-6.1", "g++-6", "g++-5.3", "g++-5.2", "g++-5.1", "g++-5", "g++-4.9", "g++-4.8", "clang++", "clang++-3.8", "clang++-3.7"]
stdflag = None

#######################################
# read the build.conf IFF it exists
#######################################
if os.path.isfile("build.conf"):
	import configparser
	confPath = "build.conf"

	conf = configparser.ConfigParser()
	conf.read([confPath])     # read the configuration file

	cppComp = conf.get("compiler", "cpp")
	if not cppComp in candidates:
		# insert specified compiler from build.conf at the beginning
		candidates.insert(0, cppComp)
	else:
		# move candidate to the beginning
		candidates.insert(0, candidates.pop(candidates.index(cppComp)))

	## C++14 support
	if stdflag is None:
		try:
			stdflag = conf.get("compiler", "std14")
		except:
			pass


#######################################
# determine and set compiler or exit if there is no suitable compiler
#######################################
# temporarily disable compiler check on windows.
if not sys.platform == 'Windows':
	# check CXX environment variable for default C++ compiler
	try:
		default_candidate = os.environ["CXX"]
		if not default_candidate in candidates:
			# insert specified compiler from build.conf at the beginning
			candidates.insert(0, default_candidate)
		else:
			# move candidate to the beginning
			candidates.insert(0, candidates.pop(candidates.index(default_candidate)))
	except:
		pass
	# check if the specified compiler is suitable
	if stdflag is None or len(stdflag) == 0:
		cppcompiler, stdflag = determineCompiler(candidates, ["c++14","c++11"])
	else:
		cppcompiler,_ = determineCompiler(candidates, [stdflag])
	if cppcompiler is not None:
		os.environ["CC"] = cppcompiler
		os.environ["CXX"] = cppcompiler
	else:
		errorMessages.append("ERROR: Test compilation with the following binary names was unsuccessful: {}. Make sure you have either g++ (>= 4.8) or clang++ (>= 3.7) properly installed".format(", ".join(candidates)))
		abort_installation = True

# early abort installation in case the compiler requirements aren't satisfied
if abort_installation:
	for msg in errorMessages:
		print(msg)
	exit(1)


################################################
# get the optional arguments for the compilation
################################################
parser = ArgumentParser()
parser.add_argument("-j", "--jobs", dest="jobs", help="specify number of jobs")
parser.add_argument("-o", "--optimize", dest="optimize", help="specify build type: Opt=optimize, Dbg=debug, Pro=profiling")
parser.add_argument("-c", "--with-cpp-tests", dest="cpptests", help="Also compile and run the C++ unit tests",action='store_true')
(options,args) = parser.parse_known_args()

# set optional arguments to parsed ones or the default ones
if options.jobs is not None:
	jobs = options.jobs
else:
	jobs = multiprocessing.cpu_count()
if options.optimize is not None:
	optimize = options.optimize
else:
	optimize = "Opt"

# make sure sys.argv is correct for setuptools
sys.argv = [__file__] + args


# this defintion probably has to stand here...
def build_NetworKit():
	if scons_available:
		comp_cmd = ["scons", "--optimize={0}".format(optimize), "--target=Core", "-j{0}".format(jobs)]
		# scons is available, now check if the user has created a build.conf
		if not os.path.isfile("build.conf"):
			# we assume, we're in a clone of the repository or in an archived copy of the repository and the user/developer has NOT created a build.conf
			# and therefore needs the information about the compiler
			comp_cmd.append("--compiler={0}".format(cppcompiler))
			comp_cmd.append("--std={0}".format(stdflag))
		print("initializing NetworKit compilation with: {0}".format(" ".join(comp_cmd)))
		if not subprocess.call(comp_cmd) == 0:
			print("scons returned an error, exiting setup.py")
			exit(1)
	else:
		from mbe import MinimalBuildEnvironment
		# minimal builder as fallback for scons
		def_compile_flags = ["-c", "-std={}".format(stdflag), "-Wall", "-fmessage-length=0", "-fPIC", "-fopenmp"]
		release_compile_flags = ["-O3", "-DNDEBUG", "-DLOG_LEVEL=LOG_LEVEL_INFO"]
		builder = MinimalBuildEnvironment(def_compile_flags,"",release_compile_flags,"","Opt", cppcompiler, "networkit/cpp")
		builder.compile("Core")

# this defintion probably has to stand here...
def additional_clean():
	clean_cmd = ["scons", "--optimize={0}".format(optimize), "--target=Core", "-c"]
	subprocess.call(clean_cmd)
	if cython_available and os.path.isfile("networkit/_NetworKit.cpp") and "clean" in sys.argv:
		os.remove("networkit/_NetworKit.cpp")


#class CustomCleanCmd(CleanCmd):
#	def initialize_options(self):
#		CleanCmd.initialize_options(self)

#	def finalize_options(self):
#		CleanCmd.finalize_options(self)

#	def run(self):
#		additional_clean()
#		CleanCmd.run(self)

class CustomStBuildExtCmd(SetuptoolsBuildExtCmd):
	def initialize_options(self):
		SetuptoolsBuildExtCmd.initialize_options(self)

	def finalize_options(self):
		SetuptoolsBuildExtCmd.finalize_options(self)

	def run(self):
		build_NetworKit()
		SetuptoolsBuildExtCmd.run(self)

class CustomCythonBuildExtCmd(CythonBuildExtCmd):
	def initialize_options(self):
		CythonBuildExtCmd.initialize_options(self)

	def finalize_options(self):
		CythonBuildExtCmd.finalize_options(self)

	def run(self):
		build_NetworKit()
		if os.path.isfile("networkit/_NetworKit.cpp"):
			os.remove("networkit/_NetworKit.cpp")
		CythonBuildExtCmd.run(self)

class MyTestCommand(TestCommand):
	def initialize_options(self):
		TestCommand.initialize_options(self)

	def finalize_options(self):
		TestCommand.finalize_options(self)

	def run(self):
		if options.cpptests:
			optimize = "Dbg"
			comp_cmd = ["scons", "--optimize={0}".format(optimize), "--target=Tests", "-j{0}".format(jobs)]
			print("initializing NetworKit compilation with: {0}".format(" ".join(comp_cmd)))
			if not subprocess.call(comp_cmd) == 0:
				print("scons returned an error, exiting setup.py")
				exit(1)
			run_cpp_cmd = ["./NetworKit-Tests-{0}".format(optimize),"-t"]
			if subprocess.call(run_cpp_cmd) == 0:
				print("C++ unit tests didn't report any errors")
			else:
				print("some C++ unit tests failed, see above")
		TestCommand.run(self)

class CustomInstallCmd(InstallCmd):
	def initialize_options(self):
		InstallCmd.initialize_options(self)

	def finalize_options(self):
		InstallCmd.finalize_options(self)

	def run(self):
		# run setuptools install command
		InstallCmd.run(self)
		# collect and print warnings about external packages used by NetworKit
		warnMessages = collectExternalPackageStatus()
		if len(warnMessages) > 0:
			for msg in warnMessages:
				print(msg)
			print("Save this list and check for each package how to install it on your system.")


src = []
# src can either be _NetworKit.pyx or the cythonized _NetworKit.cpp
do_cythonize = False
# depending on the role in which the setup script is called, it will be determined if _NetworKit.pyx will be cythonized.
build_ext_cmd = None
# the `build_ext` command depends on the role of the setup script
if not os.path.exists(".hg") and os.path.isfile("networkit/_NetworKit.cpp"):
	#print("using pre-cythonized _NetworKit.cpp")
	# we assume, were not in the repository, but installing networkit from a zip or via pip
	if os.path.isfile("MANIFEST.in"):
		os.remove("MANIFEST.in")
	src = ["networkit/_NetworKit.cpp"]
	build_ext_cmd = CustomStBuildExtCmd
elif os.path.isfile("networkit/_NetworKit.pyx") and cython_available:
	#print("cythonize _NetworKit.pyx to _NetworKit.cpp")
	# remove _NetworKit.cpp to make room for cython
	#if cython_available and os.path.isfile("networkit/_NetworKit.cpp"):
	#	os.remove("networkit/_NetworKit.cpp")
	build_ext_cmd = CustomCythonBuildExtCmd
	src = ["networkit/_NetworKit.pyx"]
	do_cythonize = True
else:
	print("ERROR: Some requirements aren't met.\nIf you try to install/build NetworKit from a clone of the repository or a ZIP archive, make sure you have Cython (version >= 0.21) installed under the __same__ Python 3 version from which you tried to install NetworKit.\nExiting...""")
	exit(1)

# initialize Extension module with the appropriate source file
modules = [Extension("_NetworKit",
	src,
	language = "c++",
	extra_compile_args=["-fopenmp", "-std={}".format(stdflag), "-O3", "-DNOGTEST"],
	extra_link_args=["-fopenmp", "-std={}".format(stdflag)],
	libraries=["NetworKit-Core-{0}".format(optimize)],
	library_dirs=["./"])]

if do_cythonize:
	for e in modules:
		e.cython_directives = {"embedsignature" : True}

# initialize the setup with the appropriate commands.
setup(
	name			= version.name,
	version			= version.version,
	author			= version.author,
	author_email	= version.author_email,
	url				= version.url,
	download_url	= version.download_url,
	description		= version.description,
	long_description= version.long_description,
	license			= version.license,
	packages		= find_packages(),
	package_data	= {'networkit.profiling': ['html/*','latex/*','description/*']},
	keywords		= version.keywords,
	platforms		= version.platforms,
	classifiers		= version.classifiers,
	cmdclass		= {'build_ext' : build_ext_cmd, 'test' : MyTestCommand, 'install' : CustomInstallCmd}, #'clean' : CustomCleanCmd,
	test_suite		= 'nose.collector',
	ext_modules		= modules,
	zip_safe		= False)
