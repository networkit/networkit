from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import os

os.environ["CC"] = "g++-4.7"
os.environ["CXX"] = "g++-4.7"


modules = [Extension("NetworKit",
					 ["NetworKit.pyx", "../src/auxiliary/RandomInteger.cpp", "../src/graph/Graph.cpp"],
					 language = "c++",
					 extra_compile_args=["-fopenmp", "-std=c++11"],
					 extra_link_args=["-fopenmp", "-std=c++11"])]

for e in modules:
	e.cython_directives = {"embedsignature" : True}

setup(name="NetworKit",
	 cmdclass={"build_ext": build_ext},
	 ext_modules=modules)

