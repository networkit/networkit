import version
import sys
if "setuptools" not in sys.modules:
	from ez_setup import use_setuptools
	# in case setuptools is not installed
	use_setuptools()

from setuptools import setup
from setuptools import Extension
from setuptools import find_packages
#from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import multiprocessing
import os
import shutil

from subprocess import Popen
import shlex

from argparse import ArgumentParser

try:
	if shutil.which("scons") is None:
		print("ERROR: Build system SCons is not installed. Please install and rerun setup.py")
		exit(1)
except:
	print("WARNING: unable to check whether build system SCons is installed")
	

# remove MANIFEST.in when networkit is not in the repository (i.e. a download pypi/zip)
if os.path.isfile("MANIFEST.in") and not os.path.exists(".hg"):
	os.remove("MANIFEST.in")


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
gcc_versions = ["","-4.9","-4.8"] #"4.7"
gcc = ""
v = 0
while gcc_version_satisfied == False and v < len(gcc_versions):
	try:
		comp_cmd = "g++{0} -o test sample.cpp -fopenmp -std=c++11".format(gcc_versions[v])
		#print(comp_cmd)
		comp_proc = Popen(shlex.split(comp_cmd))
		comp_proc.wait()
		if (comp_proc.returncode == 0):
			gcc_version_satisfied = True
			gcc = "g++{0}".format(gcc_versions[v])
			#print("your latest gcc is {0}".format(gcc))
	except:
		foo = 0
		#print("g++-{0} is not installed".format(gcc_versions[v]))
	v += 1
os.remove("sample.cpp")
if gcc_version_satisfied:
	os.remove("test")
	os.environ["CC"] = gcc
	os.environ["CXX"] = gcc
else:
	print("please install GCC 4.8 or later to be able to install NetworKit")
	exit(1)



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
#for e in sys.argv:
#	print(e)
#print("################")

def build_NetworKit():
	#os.chdir("./networkit")
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

 #print("source files: {0}".format(src))

modules = [Extension("_NetworKit",
					src,
					language = "c++",
					extra_compile_args=["-fopenmp", "-std=c++11", "-O3", "-DNOGTEST"],
					extra_link_args=["-fopenmp", "-std=c++11"],
					libraries=["NetworKit-Core-{0}".format(optimize)],
					library_dirs=["./"])]

for e in modules:
	e.cython_directives = {"embedsignature" : True}

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
#	packages = ["networkit"],
#	package_dir ={"networkit" : "./","networkit.gephi":"networkit","networkit.viztools":"networkit"},
	packages	= find_packages(),
#	packages = ["networkit","networkit.gephi","networkit.viztools"].append(modules),
#	package_dir = packages_dir,
	keywords	= version.keywords,
	platforms	= version.platforms,
	classifiers	= version.classifiers,
	cmdclass	= {"build_ext": build_ext},
	ext_modules	= modules,
	install_requires= version.install_requires,
	zip_safe	= False)#,

 #print("[Done] setup.py")
