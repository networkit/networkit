import sys
import shutil
import os

cmakeCompiler = "" #possibilty to specify a compiler
buildDirectory = "build_python"
ninja_available = False

if sys.version_info.major < 3:
	print("ERROR: NetworKit requires Python 3.")
	sys.exit(1)
if shutil.which("cmake") is None:
	print("ERROR: NetworKit compilation requires cmake.")
	sys.exit(1)

ninja_available = shutil.which("ninja") is not None
if not ninja_available or shutil.which("make") is None:
	print("ERROR: NetworKit compilation requires Make or Ninja.")
	sys.exit(1)
try:
	from setuptools import setup # to ensure setuptools is installed
except ImportError:
	print("ERROR: Setuptools is required to install networkit python module.\nInstall via pip3 install setuptools.")
	sys.exit(1)

import os
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
else:
	import multiprocessing
	jobs = multiprocessing.cpu_count()
# make sure sys.argv is correct for setuptools
sys.argv = [__file__] + args

################################################
# compiler identification
################################################

candidates = ["g++", "g++-8", "g++-7", "g++-6.1", "g++-6", "g++-5.3", "g++-5.2", "g++-5.1", "g++-5", "g++-4.9", "g++-4.8", "clang++", "clang++-3.8", "clang++-3.7"]

def determineCompiler(candidates, std, flags):
	#TODO: proper c++11 test?
	sample = open("sample.cpp", "w")
	sample.write("""
	#include <omp.h>
	#include <iostream>
	void helloWorld() {
		std::cout << "Hello world" << std::endl;
	}
	int main (int argc, char *argv[]) {
		helloWorld();
		int nthreads, tid;
		#pragma omp parallel private(nthreads, tid)
		{
			tid = omp_get_thread_num();
			std::cout << \"Hello World from thread = \" << tid << std::endl;
			if (tid == 0) {
				nthreads = omp_get_num_threads();
				std::cout << \"Number of threads = \" << nthreads << std::endl;
			}
		}
	}""")
	sample.close()
	for compiler in candidates:
		cmd = [compiler,"-o","test_build","-std={}".format(std)]
		cmd.extend(flags)
		cmd.append("sample.cpp")
		try:
			if subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0:
				os.remove("sample.cpp")
				os.remove("test_build")
				return compiler, std
		except:
			pass
	return "",""

# only check for a compiler if none is specified
if cmakeCompiler == "":
	cmakeCompiler,_ = determineCompiler(candidates, "c++11", ["-fopenmp"])
	if cmakeCompiler == "":
		print("ERROR: No suitable compiler found. Install any of these: ",candidates)
		exit(1)

################################################
# functions for cythonizing and building networkit
################################################

def cythonizeFile(filepath):
	cython_available = shutil.which("cython") is not None
	if not cython_available and not os.path.isfile(filepath.replace("pyx","cpp")):
		print("ERROR: Neither cython nor _NetworKit.cpp is provided. Build cancelled")
		exit(1)
	if not cython_available and os.path.isfile(filepath.replace("pyx","cpp")):
		print("Cython not available, but _NetworKit.cpp provided. Continue build without cythonizing")
	else:
		print("Cythonizing _NetworKit.pyx...", flush=True)
		if not os.path.isfile(filepath):
			print("_NetworKit.pyx is not available. Build cancelled.")
			exit(1)
		comp_cmd = ["cython","-3","--cplus","-t",filepath]
		if not subprocess.call(comp_cmd) == 0:
			print("cython returned an error, exiting setup.py")
			exit(1)
		print("_NetworKit.pyx cythonized", flush=True)

def buildNetworKit(withTests = False):
	# Cythonize file
	cythonizeFile("networkit/_NetworKit.pyx")
	if not os.path.isdir(buildDirectory):
		os.makedirs(buildDirectory)
	# Build cmake call
	comp_cmd = ["cmake","-DCMAKE_BUILD_TYPE=Release"]
	comp_cmd.append("-DCMAKE_CXX_COMPILER="+cmakeCompiler)
	from sysconfig import get_paths, get_config_var
	comp_cmd.append("-DNETWORKIT_PYTHON="+get_paths()['include']) #provide python.h files
	comp_cmd.append("-DNETWORKIT_PYTHON_SOABI="+get_config_var('SOABI')) #provide lib env specification
	if ninja_available:
		comp_cmd.append("-GNinja")
	comp_cmd.append("..") #call CMakeLists.txt from networkit root
	# Run cmake
	print("initializing NetworKit compilation with: '{0}'".format(" ".join(comp_cmd)), flush=True)
	if not subprocess.call(comp_cmd, cwd=buildDirectory) == 0:
		print("cmake returned an error, exiting setup.py")
		exit(1)
	build_cmd = []
	if ninja_available:
		build_cmd = ["ninja", "-j"+str(jobs)]
	else:
		build_cmd = ["make", "-j"+str(jobs)]
	print("Build with: '{0}'".format(" ".join(build_cmd)), flush=True)
	if not subprocess.call(build_cmd, cwd=buildDirectory) == 0:
		print("Build tool returned an error, exiting setup.py")
		exit(1)

################################################
# custom build commands to integrate with setuptools
################################################
from setuptools.command.test import test as TestCommand
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools import Command

class buildNetworKitCommand(_build_ext):
	def initialize_options(self):
		_build_ext.initialize_options(self)

	def finalize_options(self):
		_build_ext.finalize_options(self)

	def run(self):
		buildNetworKit()
		from sysconfig import get_config_var
		libname = "_NetworKit."+get_config_var('SOABI')+".so"
		shutil.copyfile(buildDirectory+"/"+libname,libname)
		#os.symlink(buildDirectory+"/"+libname, libname) # would be nicer but setup.py install has to be adapted then

class CustomTestCommand(TestCommand):
	def initialize_options(self):
		TestCommand.initialize_options(self)

	def finalize_options(self):
		TestCommand.finalize_options(self)

	def run(self):
		print("Testing the python module is currently not supported")
		exit(1)
		TestCommand.run(self)

################################################
# initialize python setup
################################################
from setuptools import find_packages # in addition to setup
import version
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
	cmdclass			= {'build_ext' : buildNetworKitCommand, 'test' : CustomTestCommand},
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
