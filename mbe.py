import os
import subprocess
import shlex
import multiprocessing
import fnmatch
from concurrent.futures import ProcessPoolExecutor

class MinimalBuildEnvironment:
	""" A very minimalistic Build Environment that aims to replace SCons for the
		NetworKit installation via PIP/easy_install on systems where SCons is not available. 
		Currently, it is only possible to compile NetworKit as a library in optimization flags, 
		which is the desired behavior for installation.
		It remains unclear if this build environment will be extended to fully replace SCons.
	"""

	def __init__(self, default_compile_flags,debug_compile_flags,release_compile_flags,linker_flags,optimize,compiler, src_dir):
		""" Constructor. """
		self.__default_compile_flags	= default_compile_flags
		self.__debug_compile_flags		= debug_compile_flags
		self.__release_compile_flags	= release_compile_flags
		self.__linker_flags				= linker_flags
		self.__optimize					= optimize
		self.__build_dir				= ".build{0}".format(optimize)
		self.__compiler					= compiler
		self.__src_dir					= src_dir
		self.__object_files				= []

	def __getSourceFiles(self, target):
		""" This functions gathers and returns all .cpp files for the given target. """
		source = []

		# walk source directory and find ONLY .cpp files
		for (dirpath, dirnames, filenames) in os.walk(self.__src_dir):
			for name in fnmatch.filter(filenames, "*.cpp"):
				source.append(os.path.join(dirpath, name))

		# exclude files depending on target, executables will be addes later
		xpatterns = ["*-X.cpp"]
		excluded = []

		# only the target "Test" requires Benchmark and GTest files
		if (target not in ["Tests"]):
			# exclude files matching following patterns
			xpatterns += ["*GTest.cpp","*Benchmark.cpp"]

		for pattern in xpatterns:
			for name in fnmatch.filter(source, pattern):
				excluded.append(name)

		#print("excluded source files: {0}".format(excluded))
		source = [name for name in source if name not in excluded]

		# add executable
		if target == "Tests":
			source.append(os.path.join(srcDir, "Unittests-X.cpp"))
		elif target in ["Core","Lib"]:
			pass # no executable
		else:
			print("Unknown target: {0}".format(target))
			exit(1)
		return source

	def compile_file(self, file):
		""" Compiles a given file and returns the name of the output file. """
		ofile = "{0}.o".format(file.split("/")[-1][:-4])
		comp_cmd = [self.__compiler] + self.__default_compile_flags + self.__release_compile_flags + ["-o{0}".format(os.path.join(self.__build_dir,ofile)), file]
		print(" ".join(comp_cmd))
		return (subprocess.call(comp_cmd),ofile)

	def compile(self, target):
		""" Compiles the all source files and runs the appropriate commands for the given target. """
		# make build dir if not existing yet
		if not os.path.exists(self.__build_dir):
			os.mkdir(self.__build_dir)

		# get source files
		cppfiles = self.__getSourceFiles(target)
		#print(cppfiles)
		# compile each source file on its own. halt if an error occurs
		with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
			for (returncode, ofile) in executor.map(self.compile_file,cppfiles):
				self.__object_files.append(ofile)
				if not returncode == 0:
					print("compilation of a file went wrong, exiting...")
					exit(1)

		# pull together each object file in one string for linking/archiving
		#linker_sources_str = ""
		linker_sources = []
		for o in self.__object_files:
			#linker_sources_str += (os.path.join(self.__build_dir,o) + " ")
			linker_sources.append(os.path.join(self.__build_dir,o))


		# link/archive files
		link_cmd = ["ar","rc", "libNetworKit-Core-{0}.a".format(self.__optimize)] + linker_sources
		print(" ".join(link_cmd))
		if not subprocess.call(link_cmd) == 0: #(link_proc.returncode != 0):
			print("error during linking/archiving, exiting...")
			exit(1)
		# index archive (if this build environment ever gets extended, this should be changed...)
		index_cmd = ["ranlib","libNetworKit-Core-{0}.a".format(self.__optimize)]
		print(" ".join(index_cmd))
		if not subprocess.call(index_cmd) == 0:
			print("error during indexing, exiting...")
			exit(1)

if __name__ == "__main__":
	DEFAULTCOMPILEFLAGS = ["-c", "-std=c++11", "-Wall", "-fmessage-length=0", "-fPIC", "-fopenmp"]
	DEBUGCOMPILEFLAGS = ["-O0", "-g3", "-DLOG_LEVEL=LOG_LEVEL_TRACE"]
	RELEASECOMPILEFLAGS = ["-O3", "-DNDEBUG", "-DLOG_LEVEL=LOG_LEVEL_INFO"]
	LINKERFLAGS = ""

	# test functionality by compiling NetworKit as a library.
	builder = MinimalBuildEnvironment(DEFAULTCOMPILEFLAGS,"",RELEASECOMPILEFLAGS,"","Opt","g++", "networkit/cpp")
	builder.compile("Core")