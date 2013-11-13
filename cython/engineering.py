""" 
Tools for algorithm engineering.
"""

from _NetworKit import setNumberOfThreads, getMaxNumberOfThreads, getCurrentNumberOfThreads

import stopwatch


def strongScaling(algorithm, threadSequence, input, inputTitle=None, repetitions=1):
	""" Evaluate strong scaling, i.e. how the performance varies with the number of threads
		for a fixed input size.
	"""
	data = []	# collects data about the experiments
	threadsAvailable = getMaxNumberOfThreads()	# remember maximum number of threads and restore later
	for nThreads in threadSequence:
		setNumberOfThreads(nThreads)
		print("set number of threads to {0}".format(getMaxNumberOfThreads()))
		for r in range(repetitions):
			print("running {0}".format(algorithm.toString()))
			timer = stopwatch.Timer()
			result = algorithm.run(input)
			timer.stop()
			print("elapsed time: {0}".format(timer.elapsed))
			if inputTitle is None:
				try:
					inputTitle = input.toString()
				except:
					inputTitle = str(input)
			# append run data
			data.append({"algo": algorithm.toString(), "input": inputTitle, "t": timer.elapsed})
	setNumberOfThreads(threadsAvailable)
	return data

def weakScaling(algorithm, threadSequence, inputSequence, inputTitles=None, repetitions=1):
	""" Evaluate weak scaling, i.e. how the performance varies with the number of threads
		for a fixed input size per processor.
	"""
	threadsAvailable = getMaxNumberOfThreads()	# remember maximum number of threads and restore later
	for (input, nThreads) in zip(inputSequence, threadSequence):
		setNumberOfThreads(nThreads)
		print("set number of threads to {0}".format(getMaxNumberOfThreads()))
		for r in range(repetitions):
			print("running {0}".format(algorithm.toString()))
			result = algorithm.run(input)
	setNumberOfThreads(threadsAvailable)
