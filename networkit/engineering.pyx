"""
Tools for algorithm engineering.
"""
# local imports
from . import stopwatch

# external imports
import csv
import warnings

cdef extern from "<networkit/auxiliary/Parallelism.hpp>" namespace "Aux":

	void _setNumberOfThreads "Aux::setNumberOfThreads" (int)
	int _getCurrentNumberOfThreads "Aux::getCurrentNumberOfThreads" ()
	int _getMaxNumberOfThreads "Aux::getMaxNumberOfThreads" ()
	void _enableNestedParallelism "Aux::enableNestedParallelism" ()

def setNumberOfThreads(nThreads):
	""" Set the number of OpenMP threads """
	_setNumberOfThreads(nThreads)

def getCurrentNumberOfThreads():
	""" Get the number of currently running threads"""
	return _getCurrentNumberOfThreads()

def getMaxNumberOfThreads():
	""" Get the maximum number of available threads"""
	return _getMaxNumberOfThreads()

def enableNestedParallelism():
	""" Enable nested parallelism for OpenMP"""
	from warnings import warn
	warn("Nested parallelism has been deprecated.")

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
