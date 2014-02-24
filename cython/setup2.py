from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import os
import shutil

from subprocess import Popen
import shlex
import sys

if "build_ext" in sys.argv:
	os.chdir("../")
	try:
		jobs = sys.argv[2]
	except:
		print("number of threads for compilation set to 4")
		jobs = "4"
	try:
		optimize = sys.argv[3]
	except:
		print("compilation mode set to optimize")
		optimize = "O"
	comp_cmd = "scons --optimize="+optimize+" --target=Core -j"+jobs
	print("initalizing NetworKit compilation with: "+comp_cmd)
	comp_proc = Popen(shlex.split(comp_cmd))
	comp_proc.wait()
	os.chdir("./cython")
elif "clean" in sys.argv:
	try:
		optimize = sys.argv[2]
	except:
		print("optimized NetworKit compilation will be cleaned")
		optimize = "O"
	os.chdir("../")
	clean_cmd = "scons --optimize="+optimize+" --target=Core -c"
	clean_proc = Popen(shlex.split(clean_cmd))
	clean_proc.wait()
	os.chdir("./cython")
	#os.rmdir("./build")
	os.remove("_NetworKit.cpp")
	
# make sure sys.argv is correct for distutils:
sys.argv = sys.argv[0:2]

# try-catch block when shutil.which is not available
try:
	if (shutil.which("g++-4.8") is not None):
		os.environ["CC"] = "g++-4.8"
		os.environ["CXX"] = "g++-4.8"

	elif (shutil.which("g++-4.7") is not None):
		os.environ["CC"] = "g++-4.7"
		os.environ["CXX"] = "g++-4.7"

	else:
		pass
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
					libraries=["NetworKit-Core-"+optimize],
					library_dirs=["../"])]

for e in modules:
	e.cython_directives = {"embedsignature" : True}

setup(name="NetworKit",
	author="Christian L. Staudt",
	author_email = "christian.staudt@kit.edu",
	description = "NetworKit is a toolbox for high-performance network analysis",
	license = "MIT",
	keywords = "graph algorithm network analysis social network",
	version="2.1",
	cmdclass={"build_ext": build_ext},
	ext_modules=modules,
	py_modules = ["NetworKit.py"])

print("[Done] setup.py")
