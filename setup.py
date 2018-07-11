import sys
import shutil
import os

initialPwd = os.getcwd()
cmakeCompiler = "" #possibilty to specify a compiler
#buildDirectory = os.path.join(initialPwd, "buildPython") #directory to build networkit to
buildDirectory = "buildPython"

if sys.version_info.major < 3:
	print("ERROR: NetworKit requires Python 3.")
	sys.exit(1)
if shutil.which("cmake") is None:
	print("ERROR: NetworKit compilation requires cmake.")
	sys.exit(1)
if shutil.which("make") is None:
	print("ERROR: NetworKit compilation requires make.")
	sys.exit(1)

from setuptools import setup # to ensure setuptools is installed

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
	comp_cmd.append("-DNETWORKIT_PYTHON=ON") # link python
	comp_cmd.append("..") #call CMakeLists.txt from networkit root
	# Run cmake
	print("initializing NetworKit compilation with: '{0}'".format(" ".join(comp_cmd)), flush=True)
	os.chdir(buildDirectory) #switch to build dir
	if not subprocess.call(comp_cmd) == 0:
		print("cmake returned an error, exiting setup.py")
		exit(1)
	# Run make
	make_cmd = ["make", "-j"+str(jobs)]
	print("run make with: '{0}'".format(" ".join(make_cmd)), flush=True)
	if not subprocess.call(make_cmd) == 0:
		print("Make returned an error, exiting setup.py")
		exit(1)
	os.chdir("../")
	#os.chdir(initialPwd)

################################################
# custom build commands to integrate with setuptools
################################################
from setuptools.command.test import test as TestCommand

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

buildNetworKit()
from sysconfig import get_config_var

shutil.copyfile(buildDirectory+"/libnetworkit.so","_NetworKit."+get_config_var('SOABI')+".so")


# initialize the setup with the appropriate commands.
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
	cmdclass			= {'test' : CustomTestCommand},
	test_suite			= 'nose.collector',
	install_requires	= version.install_requires,
	zip_safe			= False) # see https://cython.readthedocs.io/en/latest/src/reference/compilation.html

################################################
# check for tkinter installation
# cant be handled by install_requires because
# tkinter is not part of pip
################################################
try:
	import _tkinter
	del _tkinter
	import tkinter #python3 only
	del tkinter#python3 only
except:
	print("WARNING: _tkinter is necessary for NetworKit.\nPlease install _tkinter")
