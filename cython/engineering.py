""" 
Tools for algorithm engineering.
"""

from _NetworKit import setNumberOfThreads, getMaxNumberOfThreads, getCurrentNumberOfThreads


def strongScaling(algorithm, input, nThreads, repetitions):
	""" Evaluate strong scaling, i.e. how the performance varies with the number of threads
		for a fixed input size
	"""
	for nt in nThreads:
		setNumberOfThreads(nt)
		print("set number of threads to {0}".format(getMaxNumberOfThreads()))
		for r in range(repetitions):
			print("running {0}".format(algorithm.toString()))
			result = algorithm.run(input)

