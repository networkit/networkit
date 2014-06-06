import sys
if "setuptools" not in sys.modules:
	from ez_setup import use_setuptools
	# in case setuptools is not installed
	use_setuptools()

from setuptools import setup
from setuptools import Extension
#from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import multiprocessing
import os
import shutil

from subprocess import Popen
import shlex

from argparse import ArgumentParser

if shutil.which("scons") is None:
	print("Build system SCons is not installed. Please install and rerun setup.py")
	exit(1)



# get the optional arguments for the compilation
parser = ArgumentParser()
parser.add_argument("-j", "--jobs", dest="jobs", help="specify number of jobs")
parser.add_argument("-o", "--optimize", dest="optimize", help="specify build type: O=optimize, D=debug, P=profiling")
(options,args) = parser.parse_known_args()

# set optional arguments to parsed ones or the default ones
if options.jobs != None:
	jobs = options.jobs
else:
	jobs = multiprocessing.cpu_count()
if options.optimize != None:
	optimize = options.optimize
else:
	optimize = "O"

# make sure sys.argv is correct for setuptools
args.reverse()
args.append(__file__)
args.reverse() # this is not a very nice way to do this for sure
sys.argv = args
#for e in sys.argv:
#	print(e)
#print("################")

def build_NetworKit():
	#os.chdir("../")
	comp_cmd = "scons --optimize={0} --target=Core -j{1}".format(optimize,jobs)
	print("initializing NetworKit compilation with: {0}".format(comp_cmd))
	comp_proc = Popen(shlex.split(comp_cmd))
	comp_proc.wait()
	if (comp_proc.returncode != 0):
		print("scons returned an error, exiting setup.py")
		exit(1)
	os.chdir("./src/python")
	try:
		os.remove("_NetworKit.cpp")
	except:
		print("_NetworKit.cpp already deleted")

def additional_clean():
	#os.chdir("../")
	clean_cmd = "scons --optimize={0} --target=Core -c".format(optimize)
	clean_proc = Popen(shlex.split(clean_cmd))
	clean_proc.wait()
	os.chdir("./src/python")
	#os.rmdir("./build")
	try:
		os.remove("_NetworKit.cpp")
	except:
		print("_NetworKit.cpp already deleted")



if ("build_ext" in sys.argv):
	build_NetworKit()
elif (("develop" in sys.argv) and ("--uninstall" not in sys.argv)):
	try:
		os.mkdir("src/python/NetworKit")
	except:
		foo = 0
	build_NetworKit()
elif "clean" in sys.argv:
	additional_clean()
	
# try-catch block when shutil.which is not available
try:
	if (shutil.which("g++-4.8") is not None):
		os.environ["CC"] = "g++-4.8"
		os.environ["CXX"] = "g++-4.8"

	elif (shutil.which("g++-4.7") is not None):
		os.environ["CC"] = "g++-4.7"
		os.environ["CXX"] = "g++-4.7"

	else:
		os.environ["CC"] = "g++"
		os.environ["CXX"] = "g++"
except:
	os.environ["CC"] = "g++"
	os.environ["CXX"] = "g++"

print("Using compilers: {0} and {1}".format(os.environ["CC"], os.environ["CXX"]))

src = ["_NetworKit.pyx"]	# list of source files
			
print("source files: {0}".format(src))

modules = [Extension("_NetworKit",
					src,
					language = "c++",
					extra_compile_args=["-fopenmp", "-std=c++11", "-O3", "-DNOGTEST"],
					extra_link_args=["-fopenmp", "-std=c++11"],
					libraries=["NetworKit-Core-{0}".format(optimize)],
					library_dirs=["../../"])]

for e in modules:
	e.cython_directives = {"embedsignature" : True}

setup(name="NetworKit",
	author="Christian L. Staudt, Henning Meyerhenke",
	author_email = "christian.staudt@kit.edu, meyerhenke@kit.edu",
	description = "NetworKit is a toolbox for high-performance network analysis",
	license = "MIT",
	keywords = "graph algorithm network analysis social network",
	version="3.1",
	cmdclass={"build_ext": build_ext},
	ext_modules=modules,
	py_modules = ["NetworKit.py"])

print("[Done] setup.py")
