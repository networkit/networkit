from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import os
import shutil

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

srcDir = "../src"
src = ["_NetworKit.pyx"]	# list of source files
			
print("source files: {0}".format(src))

modules = [Extension("_NetworKit",
					src,
					language = "c++",
					extra_compile_args=["-fopenmp", "-std=c++11", "-O3", "-DNOGTEST"],
					extra_link_args=["-fopenmp", "-std=c++11"],
					libraries=["NetworKit-Core-O", "log4cxx", "tbb"],
					library_dirs=["../", "../Core-O/"])]

for e in modules:
	e.cython_directives = {"embedsignature" : True}

setup(name="_NetworKit",
	 cmdclass={"build_ext": build_ext},
	 ext_modules=modules)

print("[Done] setup.py")
