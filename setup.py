#!/usr/bin/python3

import shutil
import sys
import sysconfig
import os
from Cython.Build import cythonize
import numpy as np

cmakeCompiler = None
buildDirectory = "build/build_python"
ninja_available = False
enable_osx_crossbuild = False
build_tests = False

if sys.version_info.major < 3:
	print("ERROR: NetworKit requires Python 3.")
	sys.exit(1)

if "CXX" in os.environ:
	cmakeCompiler = os.environ["CXX"]

if "NETWORKIT_OVERRIDE_CXX" in os.environ:
	cmakeCompiler = os.environ["NETWORKIT_OVERRIDE_CXX"]

if "NETWORKIT_OSX_CROSSBUILD" in os.environ:
	enable_osx_crossbuild = True

if "NETWORKIT_BUILD_TESTS" in os.environ:
	build_tests = True

if shutil.which("cmake") is None:
	print("ERROR: NetworKit compilation requires cmake.")
	sys.exit(1)

ninja_available = shutil.which("ninja") is not None
if not (ninja_available or shutil.which("make")):
	print("ERROR: NetworKit compilation requires Make or Ninja.")
	sys.exit(1)
try:
	from setuptools import setup # to ensure setuptools is installed
except ImportError:
	print("ERROR: Setuptools is required to install networkit python module.\nInstall via pip3 install setuptools.")
	sys.exit(1)

os_soabi = sysconfig.get_config_var('SOABI')
if os_soabi is None:
	os_soabi =  sysconfig.get_config_var('EXT_SUFFIX').split(".")[1] # get_config_var('SOABI') is None on win32-systems

import os
import platform
import subprocess #calling cmake, make and cython

################################################
# get the optional arguments for the compilation
################################################
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-j", "--jobs", dest="jobs", help="specify number of jobs")
(options,args) = parser.parse_known_args()
if options.jobs is not None:
	jobs = options.jobs
if "NETWORKIT_PARALLEL_JOBS" in os.environ:
	jobs = int(os.environ["NETWORKIT_PARALLEL_JOBS"])
else:
	import multiprocessing
	jobs = multiprocessing.cpu_count()
# make sure sys.argv is correct for setuptools
sys.argv = [__file__] + args

################################################
# functions for cythonizing and building networkit
################################################

def buildNetworKit(install_prefix, externalCore=False, externalTlx=None, withTests=False, rpath=None, gpu=False):
	try:
		os.makedirs(buildDirectory)
	except FileExistsError:
		pass
	# Build cmake call
	abs_prefix = os.path.join(os.getcwd(), install_prefix)
	comp_cmd = ["cmake","-DCMAKE_BUILD_TYPE=Release"]
	comp_cmd.append("-DCMAKE_INSTALL_PREFIX="+abs_prefix)
	if cmakeCompiler:
		comp_cmd.append("-DCMAKE_CXX_COMPILER="+cmakeCompiler)
	comp_cmd.append("-DNETWORKIT_FLATINSTALL=ON")
	from sysconfig import get_paths, get_config_var
	# The following cmake parameters set Python-variables. This is done to avoid differences between the 
	# python-toolchain calling setup.py and cmake-based find-mechanisms.
	comp_cmd.append("-DNETWORKIT_PYTHON="+get_paths()['include']) # provide python.h files
	comp_cmd.append("-DNETWORKIT_PYTHON_EXECUTABLE="+sys.executable) # provide cmake with Python interpreter
	comp_cmd.append("-DNETWORKIT_PYTHON_SOABI="+os_soabi) # provide lib env specification
	comp_cmd.append("-DNETWORKIT_PYTHON_VERSION="+sysconfig.get_python_version())
	comp_cmd.append("-DNETWORKIT_NUMPY="+np.get_include()) # provide numpy.h files
	if externalCore:
		if sys.platform == "win32":
			# Reasoning: only static builds are supported and libs+dlls must reside in the same folder
			print("Builds with an external core are not supported on Windows.") 
			exit(1)
		comp_cmd.append("-DNETWORKIT_BUILD_CORE=OFF")
	if build_tests:
		comp_cmd.append("-DNETWORKIT_BUILD_TESTS=ON")
	if sys.platform == "win32":
		comp_cmd.append("-DNETWORKIT_STATIC=ON") # Windows only supports static core builds
		comp_cmd.append("-DNETWORKIT_BUILDING_STATELIB=ON") # Adds dllexport
		comp_cmd.append("-DNETWORKIT_PYTHON_VERSION="+sysconfig.get_python_version())
	if enable_osx_crossbuild:
		comp_cmd.append("-DCMAKE_OSX_ARCHITECTURES='arm64'")
	if externalTlx:
		comp_cmd.append("-DNETWORKIT_EXT_TLX="+externalTlx)
	if ninja_available:
		comp_cmd.append("-GNinja")
	comp_cmd.append(os.getcwd()) #call CMakeLists.txt from networkit root
	if rpath:
		comp_cmd.append("-DNETWORKIT_PYTHON_RPATH="+rpath)
	if gpu:
		if shutil.which('nvidia-smi') is not None:
			comp_cmd.append("-DNETWORKIT_CUDA=ON")
		else:
			print("Install nvidia driver tools to build NetworKit with GPU-support, exiting setup.py.")
			exit(1)
	# Run cmake
	print("initializing NetworKit compilation with: '{0}'".format(" ".join(comp_cmd)), flush=True)
	if not subprocess.call(comp_cmd, cwd=buildDirectory) == 0:
		print("cmake returned an error, exiting setup.py")
		exit(1)
	build_cmd = []
	if ninja_available:
		build_cmd = ["ninja", "install", "-j"+str(jobs)]
	else:
		build_cmd = ["make", "install", "-j"+str(jobs)]
	print("Build with: '{0}'".format(" ".join(build_cmd)), flush=True)
	if not subprocess.call(build_cmd, cwd=buildDirectory) == 0:
		print("Build tool returned an error, exiting setup.py")
		exit(1)

################################################
# custom build commands to integrate with setuptools
################################################

# Unfortunately, there is not simple way to invoke external commands from setuptools.
# Thus, we replace the 'build_ext' stage of setuptools completely.
# This looks a bit messy as we need to make sure we satisfy the (undocumented) internal APIs
# of setuptools (and distutils, which setuptools is based on).

# Our build_ext command compiles NetworKit using cmake and installs it into a directory
# supplied by distutils/setuptools. distutils/setuptools picks up all .so files from that
# directory when creating a Python package so that we do not have to do any
# additional installation.

from setuptools import Command
from setuptools import Extension

class build_ext(Command):
	sep_by = " (separated by '%s')" % os.pathsep
	user_options = [
		('inplace', 'i',
			"ignore build-lib and put compiled extensions into the source " +
			"directory alongside your pure Python modules"),
		('include-dirs=', 'I',
			"list of directories to search for header files" + sep_by),
		('library-dirs=', 'L',
			"directories to search for external C libraries" + sep_by),
		('networkit-external-core', None,
			"use external NetworKit core library"),
		('external-tlx=', None,
			"absolute path to external tlx library"),
		('rpath=', 'r', "additional custom rpath references"),
		('enable-gpu', 'g',
			"build with support for GPU computing (CUDA)")
	]

	def initialize_options(self):
		# TODO: While we accept --include-dirs and --library-dirs, those options are not implemented yet.
		self.build_lib = None # Output directory for libraries
		self.build_temp = None # Temporary directory
		self.inplace = False
		self.include_dirs = None
		self.library_dirs = None
		self.networkit_external_core = False
		self.external_tlx = None
		self.rpath = None
		self.enable_gpu = False

		self.extensions = None
		self.package = None

	def finalize_options(self):
		self.set_undefined_options('build',
				('build_lib', 'build_lib'),
				('build_temp', 'build_temp'))

		self.extensions = self.distribution.ext_modules
		self.package = self.distribution.ext_package

	def get_source_files(self):
		sources = [ ]
		for subdir, dirs, files in os.walk('networkit'):
			for filename in files:
				sources.append(os.path.join(subdir, filename))
		return sources

	# Returns the full Python module name of an extension.
	# Prepends the ext_package setup() option to an extension (see distutils).
	def get_ext_fullname(self, extname):
		if self.package is not None:
			return self.package + '.' + extname
		return extname

	# Returns the file name of the DSO implementing a module (see distutils).
	def get_ext_filename(self, fullname):
		return fullname + '.' + os_soabi + '.so'

	def run(self):
		# A generic build_ext command for cmake would iterate over all self.extensions.
		# However, we know that we only want to build a single extension, i.e. NetworKit.
		prefix = self.build_lib
		if self.inplace:
			# The --inplace implementation is less sophisticated than in distutils,
			# but it should be sufficient for NetworKit.
			prefix = self.distribution.src_root or os.getcwd()
		buildNetworKit(prefix, externalCore=self.networkit_external_core, externalTlx=self.external_tlx, rpath=self.rpath, gpu=self.enable_gpu)

	def get_ext_fullpath(self, ext_name):
		"""Returns the path of the filename for a given extension.
		The file is located in `build_lib` or directly in the package
		(inplace option).
		"""
		fullname = self.get_ext_fullname(ext_name)
		modpath = fullname.split('.')
		filename = self.get_ext_filename(modpath[-1])

		if not self.inplace:
			# no further work needed
			# returning :
			#   build_dir/package/path/filename
			filename = os.path.join(*modpath[:-1]+[filename])
			return os.path.join(self.build_lib, filename)

		# the inplace option requires to find the package directory
		# using the build_py command for that
		package = '.'.join(modpath[0:-1])
		build_py = self.get_finalized_command('build_py')
		package_dir = os.path.abspath(build_py.get_package_dir(package))

		# returning
		#   package_dir/filename
		return os.path.join(package_dir, filename)

	def get_outputs(self):
		# And build the list of output (built) filenames.  Note that this
		# ignores the 'inplace' flag, and assumes everything goes in the
		# "build" tree.
		outputs = []
		for ext in self.extensions:
			outputs.append(self.get_ext_fullpath(ext.name))
		return outputs

################################################
# initialize python setup
################################################
from setuptools import find_packages # in addition to setup
import version

compiler_directives = {}
if build_tests:
	compiler_directives['linetrace'] = True
compiler_directives['language_level'] = 3

setup(
	name				= version.name,
	version				= version.version,
	author				= version.author,
	author_email		= version.author_email,
	url					= version.url,
	download_url		= version.download_url,
	description			= version.description,
	long_description	= version.long_description,
	license				= version.license,
	packages			= find_packages(),
	package_data		= {'networkit.profiling': ['html/*','latex/*','description/*']},
	keywords			= version.keywords,
	platforms			= version.platforms,
	classifiers			= version.classifiers,
	cmdclass			= {'build_ext': build_ext},
	ext_modules			= cythonize(["networkit/*pyx"], compiler_directives=compiler_directives),
	test_suite			= 'nose.collector',
	install_requires	= version.install_requires,
	zip_safe			= False) # see https://cython.readthedocs.io/en/latest/src/reference/compilation.html

################################################
# check for tkinter installation
# cant be handled by install_requires because
# tkinter is not part of pip
################################################
try:
	import tkinter #python3 only
	del tkinter #python3 only
except:
	print("WARNING: _tkinter is necessary for networkit.profiling module if not used via jupyter.\nInstall tkinter if needed")
