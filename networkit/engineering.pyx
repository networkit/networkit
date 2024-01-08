# distutils: language=c++

"""
Tools for algorithm engineering.
"""
from libc.stdint cimport uint64_t
from libcpp.string cimport string
from libcpp cimport bool as bool_t

from .structures cimport index, node

# local imports
from .helpers import stdstring, pystring

# external imports
import csv
import timeit
import warnings

cdef extern from "<networkit/Globals.hpp>" namespace "NetworKit":

	index _none "NetworKit::none"

none = _none

cdef extern from "<networkit/auxiliary/Parallelism.hpp>" namespace "Aux":

	void _setNumberOfThreads "Aux::setNumberOfThreads" (int)
	int _getCurrentNumberOfThreads "Aux::getCurrentNumberOfThreads" ()
	int _getMaxNumberOfThreads "Aux::getMaxNumberOfThreads" ()

def setNumberOfThreads(nThreads):
	""" 
	setNumberOfThreads(nThreads)

	Set the number of OpenMP threads

	Parameters
	----------
	nThreads : int
		Number of threads.

	"""
	_setNumberOfThreads(nThreads)

def getCurrentNumberOfThreads():
	""" 
	getCurrentNumberOfThreads()

	Get the number of currently running threads.
	
	Returns
	-------
	int
		Number of threads.
	"""
	return _getCurrentNumberOfThreads()

def getMaxNumberOfThreads():
	""" 
	getMaxNumberOfThreads()

	Get the maximum number of available threads
	
	
	Returns
	-------
	int
		Max number of threads.
	"""
	return _getMaxNumberOfThreads()

cdef extern from "<networkit/auxiliary/Log.hpp>" namespace "Aux":

	#void _configureLogging "Aux::configureLogging" (string loglevel)
	string _getLogLevel "Aux::Log::getLogLevel" () except +
	void _setLogLevel "Aux::Log::setLogLevel" (string loglevel) except +

def getLogLevel():
	""" 
	getLogLevel()
	
	Get the current log level.
	
	Returns
	-------
	logLevel
		The current loglevel.
	"""
	return pystring(_getLogLevel())

def setLogLevel(loglevel):
	""" 
	setLogLevel(loglevel)

	Set the current loglevel
	
	Parameters
	----------
	loglevel : str
		The new loglevel. Possible values: TRACE, DEBUG, INFO, WARN, ERROR, FATAL, QUIET
	"""
	_setLogLevel(stdstring(loglevel))

cdef extern from "<networkit/GlobalState.hpp>" namespace "NetworKit":

	void _setPrintLocation "NetworKit::GlobalState::setPrintLocation" (bool_t) except +

def setPrintLocation(flag):
	""" 
	setPrintLocation(flag)
	
	Switch locations in log statements on or off
	
	Parameters
	----------
	flag : bool
		Sets whether to also log file, function and line of code. Default: False.
	"""
	_setPrintLocation(flag)

cdef extern from "<networkit/auxiliary/Random.hpp>" namespace "Aux::Random":

	void _setSeed "Aux::Random::setSeed" (uint64_t, bool_t)

def setSeed(uint64_t seed, bool_t useThreadId):
	""" 
	setSeed(seed, useThreadId)
	
	Set the random seed that is used in NetworKit.

	Note that there is a separate random number generator per thread.

	Parameters
	----------
	seed : int
		The seed
	useThreadId : bool
		If the thread id shall be added to the seed
	"""
	_setSeed(seed, useThreadId)

def strongScaling(algorithmClass, threadSequence, inargs, inputTitle=None, repetitions=1, outPath=None):
	"""
	strongScaling(algorithmClass, threadSequence, inargs, inputTitle=None, repetitions=1, outPath=None)
	
	Evaluate strong scaling, i.e. how the performance varies with the number of threads for a fixed input size.

	Note
	----
	Algorithm is executed by calling :code:`algorithmClass(**inargs)`. See parameter for more details.

	Parameters
	----------
	algorithmClass :
		Algorithm, which should be tested.
	threadSequence : list(int)
		A list of thread numbers to run the :code:`algorithmClass` with.
	inargs : **kwargs
		Input arguments for algorithm. 
	inputTitle : str, optional
		Set a title for the output. Default: None
	repetitions : int, optional
		Number of repetitions. Default: 1
	outPath : str, optional
		File for writing the output to. Default: None
	"""
	data = []	# collects data about the experiments
	threadsAvailable = getMaxNumberOfThreads()	# remember maximum number of threads and restore later
	for nThreads in threadSequence:
		setNumberOfThreads(nThreads)
		print("set number of threads to {0}".format(getMaxNumberOfThreads()))
		for r in range(repetitions):
			algorithm = algorithmClass(**inargs)
			print("running algorithm")
			timer = timeit.default_timer()
			result = algorithm.run()
			timerElapsed = timeit.default_timer() - timer
			print("elapsed time: {0}".format(timerElapsed))
			if inputTitle is None:
				inputTitle = str(input)
			# append run data
			data.append({"input": inputTitle, "threads": nThreads, "time": timerElapsed})
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
	"""
	weakScaling(algorithmClass, inargs, threadSequence, inputSequence, inputTitles=None, repetitions=1, outPath=None)
	
	Evaluate weak scaling, i.e. how the performance varies with the number of threads for a fixed input size per processor.

	Note
	----
	Algorithm is executed by calling :code:`algorithmClass(input, **inargs)`. See parameter for more details.

	Parameters
	----------
	algorithmClass :
		Algorithm, which should be tested.
	inargs : **kwargs
		Input arguments for algorithm.
	threadSequence : list(int)
		A list of thread numbers to run the :code:`algorithmClass` with.
	inputSequence : list(networkit.Graph)
		A list of graphs. The input algorithm is evaluated against every list-member.
	inputTitles : str, optional
		Set a title for the output. Default: None
	repetitions : int, optional
		Number of repetitions. Default: 1
	outPath : str, optional
		File for writing the output to. Default: None
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
			print("running algorithm")
			timer = timeit.default_timer()
			result = algorithm.run()
			timerElapsed = timeit.default_timer() - timer
			# append run data
			data.append({"input": inputTitles[i], "threads": nThreads, "time": timerElapsed})
	setNumberOfThreads(threadsAvailable)
	if outPath:
		with open(outPath, "w") as outFile:
			columns = ["algo", "input", "threads", "time"]
			writer = csv.DictWriter(outFile, fieldnames=columns, delimiter="\t")
			writer.writeheader()
			for row in data:
				writer.writerow(row)
	return data
