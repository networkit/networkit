from NetworKit import *
from evaluation import *
from multiprocessing import Pool
from concurrent.futures import *
import time
import copy
from txtWriter import TxtResultWriter
from sqliteWriter import SqliteResultWriter
from paramAdjusting import *
from bbalgorithms import *
import string
import pickle
from types import *

# -----------------------------------------------------------------------
# The purpose of the following script is to automatically apply a set
# of backbone algorithms to a set of input graphs and to determine

# certain graph properties for evaluation.
# -----------------------------------------------------------------------

def userInputYesNo(q):
	while True:
		inp = input('' + q + ' [y/n]').lower()
		if inp == 'y':
			return True
		elif inp == 'n':
			return False
		else:
			print('Respond with [y/n].')

# Returns a set of backbone tasks that contains all combinations of algorithms and target edge percentages
# for the given graph. Argument represented as a tuple to be able to use it in executor.map()
def getBackboneTasks(args):
	#TODO: actually use those interval limits given in the input strings...
	#Unpack arguments
	graphFile, graphFormat, inputAlgorithms, targetEdgePercentages, outputDir = args

	#Load graph
	G = readGraph(graphFile, graphFormat)
	G.indexEdges()

	print("getBackboneTasks")

	#Generate a set of algorithms to apply to that graph
	algorithms = []
	for inputAlgorithm in inputAlgorithms:
		for targetEp in targetEdgePercentages:
			if inputAlgorithm.parameterType() == "FloatType":
				outputAlgorithm = parameterizeUsingDouble(G, inputAlgorithm, targetEp, 0.0, 1.0, 20)
			elif inputAlgorithm.parameterType() == "IntType":
				outputAlgorithm = parameterizeUsingInteger(G, inputAlgorithm, targetEp, 0, 10)
			else:
				outputAlgorithm = None
			name = inputAlgorithm.getShortName(targetEp)
			algorithms.append(BackboneAlgorithm(outputAlgorithm, name))

	bbOriginal = bb_Original()
	algorithms.append(BackboneAlgorithm(bbOriginal.getAlgorithmExpr(), bbOriginal.getShortName()))
	return BackboneTask(graphFile, graphFormat, algorithms, outputDir)

def main():
	# Create output dir
	outputDir = "./output/"
	if not os.path.exists(os.path.dirname(outputDir)):
		os.makedirs(os.path.dirname(outputDir))

	#Target edge percentages:
	targetEdgePercentages = [0.01, 0.02, 0.1, 0.5, 0.75, 0.9]

	#The algorithms to evaluate. See bbalgorithms
	inputAlgorithms = [
	#	bb_SimmelianBackboneNonParametric(),
	#	bb_SimmelianBackboneParametric(),
	#	bb_LocalSimilarityBackbone(),
		bb_SimmelianMultiscaleBackbone()
	]

	#The graphs to calculate backbones for.
	graphs = [
	#	("./input/Caltech36.graphml", Format.GraphML),
		("./input/kitEmail.graphml", Format.GraphML)
	#	("./input/Yale4.graphml", Format.GraphML),
	#	("./input/Virginia63.graphml", Format.GraphML),
	#	("./input/Tennessee95.graphml", Format.GraphML),
	#	("./input/us-aviation-t100-2013.graphml", Format.GraphML),
	#	("./input/LFR-1000.graph", Format.METIS),
	#	("./input/LFR-10000.graph", Format.METIS),
	#	("./input/bter-graph.dat", Format.EdgeListCommaOne)
		]

	#Get parameterized tasks (either by calculation or use existing file)
	pickleFile = "tasks.dat"
	tasks = []
	if not os.path.isfile('./' + pickleFile) or not userInputYesNo('Tasks.dat found. Use existing tasks?'):
		print("Parameterizing algorithms...")

		argzz = [(graphFile, graphFormat, inputAlgorithms, targetEdgePercentages, outputDir) for (graphFile, graphFormat) in graphs]
		for argz in argzz:
			tasks.append(getBackboneTasks(argz))

		with open(pickleFile, 'wb') as pf:
			pickle.dump(tasks, pf)

		print("Finished parameterizing algorithms. Wrote pickle file ", pickleFile)

	#Unpickle the tasks
	with open(pickleFile, 'rb') as pf:
		tasks = pickle.load(pf)

	#Apply the algorithms to the graphs
#	executor = ProcessPoolExecutor(max_workers=4)
#	taskResults = executor.map(executeTask, tasks)
	taskResults = []
	for task in tasks:
		taskResults.append(executeTask(task))

	#Output files...
	writers = [TxtResultWriter("./output"), SqliteResultWriter("./output/backbones.db")]
	for taskResult in taskResults:
		for writer in writers:
			writer.receiveResult(taskResult)

if __name__ == "__main__":
	sys.exit(main())
