import version
import sys
if "setuptools" not in sys.modules:
	from ez_setup import use_setuptools
	# in case setuptools is not installed
	use_setuptools()

from setuptools import setup
from setuptools import Extension
from setuptools import find_packages
from setuptools.command.test import test as TestCommand
import unittest
#from distutils.extension import Extension

abortInstallation = False
errorMessages = []
warnMessages = []
try:
	from Cython.Build import cythonize
	from Cython.Distutils import build_ext
except:
	abortInstallation = True
	errorMessages.append("ERROR: Cython not installed. Please install Cython and rerun")

import multiprocessing
import os
import shutil

from subprocess import Popen, DEVNULL
import shlex

from argparse import ArgumentParser

try:
	if shutil.which("scons") is None:
		errorMessages.append("ERROR: Build system SCons is not installed. Please install and rerun")
		abortInstallation = True
except:
	print("WARNING: unable to check whether build system SCons is installed")

#prepare sample.cpp file necessary to determine gcc
sample = open("sample.cpp", "w")
sample.write("""/*****************************************************************************
* File: sample.cpp
* DESCRIPTION:
*   OpenMP Example - Hello World - C/C++ Version
*   In this simple example, the master thread forks a parallel region.
*   All threads in the team obtain their unique thread number and print it.
*   The master thread only prints the total number of threads.  Two OpenMP
*   library routines are used to obtain the number of threads and each
*   thread's number.
* AUTHOR: Blaise Barney  5/99
* LAST REVISED: 04/06/05
******************************************************************************/
#include <omp.h>
#include <iostream>
int main (int argc, char *argv[]) {
	int nthreads, tid;
	/* Fork a team of threads giving them their own copies of variables */
	#pragma omp parallel private(nthreads, tid)
	{
		/* Obtain thread number */
		tid = omp_get_thread_num();
		std::cout << \"Hello World from thread = \" << tid << std::endl;
		/* Only master thread does this */
		if (tid == 0) {
			nthreads = omp_get_num_threads();
			std::cout << \"Number of threads = \" << nthreads << std::endl;
		}
	}  /* All threads join master thread and disband */
}""")
sample.close()


gcc_version_satisfied = False
gcc_versions = ["","-4.9","-4.8","-4.7"]
gcc = ""
v = 0
while gcc_version_satisfied == False and v < len(gcc_versions):
	try:
		comp_cmd = "g++{0} -o test_build sample.cpp -fopenmp -std=c++11".format(gcc_versions[v])
		#print(comp_cmd)
		comp_proc = Popen(shlex.split(comp_cmd), stdout=DEVNULL, stderr=DEVNULL)
		comp_proc.wait()
		if (comp_proc.returncode == 0):
			gcc_version_satisfied = True
			gcc = "g++{0}".format(gcc_versions[v])
			#print("your latest gcc is {0}".format(gcc))
	except:
		foo = 0
		#print("g++{0} is not installed".format(gcc_versions[v]))
	v += 1
os.remove("sample.cpp")
if gcc_version_satisfied:
	os.remove("test_build")
	os.environ["CC"] = gcc
	os.environ["CXX"] = gcc
else:
	errorMessages.append("ERROR: Please install GCC/g++ 4.8 or later and rerun")
	abortInstallation = True


# abort installation in case either Cython, Scons or the compiler requirements aren't satisfied
if abortInstallation:
	for msg in errorMessages:
		print(msg)
	exit(1)

# check for external packages and collect warning messages
try:
	import scipy
	del scipy
except:
	warnMessages.append("WARNING: SciPy is not installed; to use all of networkit, please install SciPy")
try:
	import numpy
	del numpy
except:
	warnMessages.append("WARNING: numpy is not installed; to use all of networkit, please install numpy")

try:
	import readline
	del readline
except:
	warnMessages.append("WARNING: readline is not installed; to use all of networkit, please install readline")

try:
	import matplotlib
	del matplotlib
except:
	warnMessages.append("WARNING: matplotlib is not installed; to use all of networkit, please install matplotlib")

try:
	import networkx
	del networkx
except:
	warnMessages.append("WARNING: networkx is not installed; to use all of networkit, please install networkx")

try:
	import tabulate
	del tabulate
except:
	warnMessages.append("WARNING: tabulate is not installed; to use all of networkit, please install tabulate")

# remove MANIFEST.in when networkit is not in the repository (i.e. a download pypi/zip)
if os.path.isfile("MANIFEST.in") and not os.path.exists(".hg"):
	os.remove("MANIFEST.in")

# remove _NetworKit.cpp, since it is very unlikely there is a scenario, where it's necessary to keep the file.
if os.path.isfile("networkit/_NetworKit.cpp"):
	os.remove("networkit/_NetworKit.cpp")

# get the optional arguments for the compilation
parser = ArgumentParser()
parser.add_argument("-j", "--jobs", dest="jobs", help="specify number of jobs")
parser.add_argument("-o", "--optimize", dest="optimize", help="specify build type: Opt=optimize, Dbg=debug, Pro=profiling")
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
args.reverse()
args.append(__file__)
args.reverse() # this is not a very nice way to do this for sure
sys.argv = args

def build_NetworKit():
	#os.chdir("./networkit")
	if os.path.isfile("build.conf"):
		comp_cmd = "scons --optimize={0} --target=Core -j{1}".format(optimize,jobs)
	else:
		comp_cmd = "scons --optimize={0} --target=Core --compiler={1} -j{2}".format(optimize,gcc,jobs)
	print("initializing NetworKit compilation with: {0}".format(comp_cmd))
	comp_proc = Popen(shlex.split(comp_cmd))
	comp_proc.wait()
	if (comp_proc.returncode != 0):
		print("scons returned an error, exiting setup.py")
		exit(1)
	#os.chdir("./src/python")
	if os.path.isfile("networkit/_NetworKit.cpp"):
		os.remove("networkit/_NetworKit.cpp")
	#os.chdir("../")

def additional_clean():
	#os.chdir("./networkit")
	clean_cmd = "scons --optimize={0} --target=Core -c".format(optimize)
	clean_proc = Popen(shlex.split(clean_cmd))
	clean_proc.wait()
	#os.chdir("./src/python")
	#os.rmdir("./build")
	if os.path.isfile("networkit/_NetworKit.cpp"):
		os.remove("networkit/_NetworKit.cpp")
	#os.chdir("../")



if "build_ext" in sys.argv:
	build_NetworKit()
elif ("develop" in sys.argv) and ("--uninstall" not in sys.argv):
	#try:
	#	os.mkdir("src/python/NetworKit")
	#except:
	#	foo = 0
	build_NetworKit()
elif "install" in sys.argv:
	build_NetworKit()
elif "clean" in sys.argv:
	additional_clean()

# try-catch block when shutil.which is not available
 #try:
	#if shutil.which("g++-4.9") is not None:
	#	os.environ["CC"] = "g++-4.9"
	#	os.environ["CXX"] = "g++-4.9"
	#elif shutil.which("g++-4.8") is not None:
	#	os.environ["CC"] = "g++-4.8"
	#	os.environ["CXX"] = "g++-4.8"

	#elif shutil.which("g++-4.7") is not None:
	#	os.environ["CC"] = "g++-4.7"
	#	os.environ["CXX"] = "g++-4.7"

	#else:
	#	pass
 #except:
	#os.environ["CC"] = "g++"
	#os.environ["CXX"] = "g++"
#print("Using compilers: {0} and {1}".format(os.environ["CC"], os.environ["CXX"]))

src = ["networkit/_NetworKit.pyx"]	# list of source files
modules = [Extension("_NetworKit",
	src,
	language = "c++",
	extra_compile_args=["-fopenmp", "-std=c++11", "-O3", "-DNOGTEST"],
	extra_link_args=["-fopenmp", "-std=c++11"],
	libraries=["NetworKit-Core-{0}".format(optimize)],
	library_dirs=["./"])]

for e in modules:
	e.cython_directives = {"embedsignature" : True}

class MyTestCommand(TestCommand):
	def initialize_options(self):
		#pass
		TestCommand.initialize_options(self)
		#loader = unittest.TestLoader()
		#TestCommand.test_suite = loader.discover('test')
		#TestCommand.loader = loader
	
	def finalize_options(self):
		#pass
		TestCommand.finalize_options(self)
	
	def run(self):
		jobs = multiprocessing.cpu_count()
		#if options.optimize is not None:
		#	optimize = options.optimize
		#else:
		optimize = "Dbg"
		comp_cmd = "scons --optimize={0} --target=Tests -j{1}".format(optimize,jobs)
		print("initializing NetworKit compilation with: {0}".format(comp_cmd))
		#comp_proc = Popen(shlex.split(comp_cmd), stdout=DEVNULL, stderr=DEVNULL)
		comp_proc = Popen(shlex.split(comp_cmd))
		comp_proc.wait()
		if (comp_proc.returncode != 0):
			print("scons returned an error, exiting setup.py")
			exit(1)
		run_cpp_cmd = "./NetworKit-Tests-{0} -t".format(optimize)
		#run_cpp_proc = Popen(shlex.split(run_cpp_cmd), stdout=DEVNULL, stderr=DEVNULL)
		run_cpp_proc = Popen(shlex.split(run_cpp_cmd))
		run_cpp_proc.wait()
		if run_cpp_proc.returncode == 0:
			print("C++ unit tests didn't report any errors")
		else:
			print("some C++ unit tests failed, see above")
#		print("return code from tests: {0}".format(run_cpp_proc.returncode))
#		build_ext.run(self)	
		comp_cmd = "scons --optimize=Opt --target=Core -j{1}".format(optimize,jobs)
		print("initializing NetworKit compilation with: {0}".format(comp_cmd))
		comp_proc = Popen(shlex.split(comp_cmd))
		comp_proc.wait()
		if (comp_proc.returncode != 0):
			print("scons returned an error, exiting setup.py")
			exit(1)

		TestCommand.run(self)


setup(
	name		= version.name,
	version		= version.version,
	author		= version.author,
	author_email	= version.author_email,
	url		= version.url,
	download_url	= version.download_url,
	description	= version.description,
	long_description= version.long_description,
	license		= version.license,
	packages	= find_packages(),
	keywords	= version.keywords,
	platforms	= version.platforms,
	classifiers	= version.classifiers,
	cmdclass	= {'build_ext' : build_ext, 'test' : MyTestCommand},
	test_suite	= 'nose.collector',
	ext_modules	= modules,
#	install_requires= version.install_requires,
	zip_safe	= False)

# print warnings about missing packages
if len(warnMessages) > 0 and "install" in sys.argv:
	for msg in warnMessages:
		print(msg)
	print("Save this list and check for each package how to install it on your system.")

 #print("[Done] setup.py")
