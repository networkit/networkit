import os
import sys
import subprocess
from subprocess import DEVNULL
import pip


def installExternalPythonPackages():
	os.system('pip3 install -r requirements.txt')

def checkDependencies():
	"""
	Checks if external dependecies are installed and warns the user if they are not. This is only meant for dependencies,
	the installation of which requires user input (e.g. must be sudo).
	"""
	necessaryDependencies = []
	try:
		import _tkinter
		del _tkinter
	except:
		necessaryDependencies.append("ERROR: _tkinter is necessary for NetworKit.\n"
			  "Please install _tkinter \n"
			  "Root privileges are necessary for this. \n"
			  "If you have these, the installation command is: sudo apt-get install python3-tk")
	return necessaryDependencies

def determineCompiler(candidates, stdFlags):
	""" This function tests a list of candidates, whether they are sufficient to the requirements of 
		NetworKit and focuses on C++11 and OpenMP support."""
	#prepare sample.cpp file necessary to determine gcc
	#TODO: proper c++11 test?
	#TODO: generalize variable names to "compiler" instead of "gcc"...
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

	[[deprecated("use the function body directly instead of wrapping it in a function.")]]
	void helloWorld() {
		std::cout << "Hello world" << std::endl;
	}

	int main (int argc, char *argv[]) {
		helloWorld();
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

	compiler_version_satisfied = False
	compiler = None
	stdflag = None
	v = 0
	i = 0
	while not compiler_version_satisfied and i < len(stdFlags):
		while not compiler_version_satisfied and v < len(candidates):
			#print("testing\t{}".format(candidates[v]))
			try:
				if subprocess.call([candidates[v],"-o","test_build","-std={}".format(stdFlags[i]),"-fopenmp","sample.cpp"],stdout=DEVNULL,stderr=DEVNULL) == 0:
					compiler_version_satisfied = True
					compiler = candidates[v]
					stdflag = stdFlags[i]
					#print("using {0} as C++ compiler with the {1} STD flag".format(candidates[v],stdFlags[i]))
			except:
				#print("{0} is not installed".format(candidates[v]))
				pass
			v += 1
		i += 1
		v = 0

	os.remove("sample.cpp")
	if compiler_version_satisfied:
		os.remove("test_build")
	return compiler, stdflag

