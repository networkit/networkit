# distutils: language=c++

"""
Tools for algorithm engineering.
"""
from libc.stdint cimport uint64_t
from libcpp.string cimport string
from libcpp cimport bool as bool_t

ctypedef uint64_t index
ctypedef index node

# local imports
from . import stopwatch

# external imports
import csv
import warnings

def pystring(stdstring):
	""" convert a std::string (= python byte string) to a normal Python string"""
	return stdstring.decode("utf-8")

def stdstring(pystring):
	""" convert a Python string to a bytes object which is automatically coerced to std::string"""
	pybytes = pystring.encode("utf-8")
	return pybytes

cdef extern from "<networkit/Globals.hpp>" namespace "NetworKit":

	index _none "NetworKit::none"

none = _none

cdef extern from "<networkit/auxiliary/Parallelism.hpp>" namespace "Aux":

	void _setNumberOfThreads "Aux::setNumberOfThreads" (int)
	int _getCurrentNumberOfThreads "Aux::getCurrentNumberOfThreads" ()
	int _getMaxNumberOfThreads "Aux::getMaxNumberOfThreads" ()

def setNumberOfThreads(nThreads):
	""" Set the number of OpenMP threads """
	_setNumberOfThreads(nThreads)

def getCurrentNumberOfThreads():
	""" Get the number of currently running threads"""
	return _getCurrentNumberOfThreads()

def getMaxNumberOfThreads():
	""" Get the maximum number of available threads"""
	return _getMaxNumberOfThreads()

cdef extern from "<networkit/auxiliary/Log.hpp>" namespace "Aux":

	#void _configureLogging "Aux::configureLogging" (string loglevel)
	string _getLogLevel "Aux::Log::getLogLevel" () except +
	void _setLogLevel "Aux::Log::setLogLevel" (string loglevel) except +
	void _setPrintLocation "Aux::Log::Settings::setPrintLocation" (bool_t) except +

def getLogLevel():
	""" Get the current log level"""
	return pystring(_getLogLevel())

def setLogLevel(loglevel):
	""" Set the current loglevel"""
	_setLogLevel(stdstring(loglevel))

def setPrintLocation(flag):
	""" Switch locations in log statements on or off"""
	_setPrintLocation(flag)

cdef extern from "<networkit/auxiliary/Random.hpp>" namespace "Aux::Random":

	void _setSeed "Aux::Random::setSeed" (uint64_t, bool_t)

def setSeed(uint64_t seed, bool_t useThreadId):
	""" Set the random seed that is used in NetworKit.

	Note that there is a separate random number generator per thread.

	Parameters:
	-----------
	seed : uint64_t
		The seed
	useThreadId : bool
		If the thread id shall be added to the seed
	"""
	_setSeed(seed, useThreadId)

def strongScaling(algorithmClass, threadSequence, inargs, inputTitle=None, repetitions=1, outPath=None):
	""" Evaluate strong scaling, i.e. how the performance varies with the number of threads
		for a fixed input size.
	"""
	data = []	# collects data about the experiments
	threadsAvailable = getMaxNumberOfThreads()	# remember maximum number of threads and restore later
	for nThreads in threadSequence:
		setNumberOfThreads(nThreads)
		print("set number of threads to {0}".format(getMaxNumberOfThreads()))
		for r in range(repetitions):
			algorithm = algorithmClass(**inargs)
			print("running {0}".format(algorithm.toString()))
			timer = stopwatch.Timer()
			result = algorithm.run()
			timer.stop()
			print("elapsed time: {0}".format(timer.elapsed))
			if inputTitle is None:
				try:
					inputTitle = input.toString()
				except:
					inputTitle = str(input)
			# append run data
			data.append({"algo": algorithm.toString(), "input": inputTitle, "threads": nThreads, "time": timer.elapsed})
	setNumberOfThreads(threadsAvailable)
	if outPath:
		with open(outPath, "w") as outFile:
			columns = ["algo", "input", "threads", "time"]
			writer = csv.DictWriter(outFile, fieldnames=columns, delimiter="\t")
			writer.writeheader()
			for row in data:
				writer.writerow(row)
	return data

def weakScaling(algorithmClass, inargs, threadSequence, inputSequence, inputTitles=None, repetitions=1, outPath=None):
	""" Evaluate weak scaling, i.e. how the performance varies with the number of threads
		for a fixed input size per processor.
	"""
	data = []	# collects data about the experiments
	threadsAvailable = getMaxNumberOfThreads()	# remember maximum number of threads and restore later
	i = -1
	for (input, nThreads) in zip(inputSequence, threadSequence):
		i += 1
		setNumberOfThreads(nThreads)
		print("set number of threads to {0}".format(getMaxNumberOfThreads()))
		for r in range(repetitions):
			algorithm = algorithmClass(input, **inargs)
			print("running {0}".format(algorithm.toString()))
			timer = stopwatch.Timer()
			result = algorithm.run()
			timer.stop()
			# append run data
			data.append({"algo": algorithm.toString(), "input": inputTitles[i], "threads": nThreads, "time": timer.elapsed})
	setNumberOfThreads(threadsAvailable)
	if outPath:
		with open(outPath, "w") as outFile:
			columns = ["algo", "input", "threads", "time"]
			writer = csv.DictWriter(outFile, fieldnames=columns, delimiter="\t")
			writer.writeheader()
			for row in data:
				writer.writerow(row)
	return data
