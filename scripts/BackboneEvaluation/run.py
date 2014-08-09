from NetworKit import *
from evaluation import *
from multiprocessing import Pool
from concurrent.futures import *
import time
import copy
from txtWriter import TxtResultWriter
from sqliteWriter import SqliteResultWriter
import string
import pickle
import parameters
import parameterization
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

#Get parameterized tasks (either by calculation or use existing file)
def getTasks(targetEdgeRatios, targetAlgorithms, targetGraphs, outputDir):
	pickleFile = "_tasks.txt"
	tasks = []
	if not os.path.isfile('./' + pickleFile) or not userInputYesNo(pickleFile + ' found. Use existing tasks?'):
		print("Parameterizing algorithms...")
		for targetGraph in targetGraphs:
			task = parameterization.getBackboneTask(targetGraph, targetAlgorithms, targetEdgeRatios, outputDir)
			tasks.append(task)

		with open(pickleFile, 'wb') as pf:
			pickle.dump(tasks, pf)

		print("Finished parameterizing algorithms. Wrote pickle file ", pickleFile)

	#Unpickle the tasks
	with open(pickleFile, 'rb') as pf:
		tasks = pickle.load(pf)

	return tasks

def writeResults(taskResults):
	#TxtResultWriter("./output"),
	writers = [SqliteResultWriter("./output/backbones.db")]
	for taskResult in taskResults:
		for writer in writers:
			writer.receiveResult(taskResult)

def main():
	# Create output dir
	outputDir = "./output/"
	if not os.path.exists(os.path.dirname(outputDir)):
		os.makedirs(os.path.dirname(outputDir))

	#Get input parameters
	edgeRatios = parameters.getEdgeRatios()
	algorithms = parameters.getAlgorithms()
	graphs = parameters.getGraphs()

	#Generate tasks from input parameters
	tasks = getTasks(edgeRatios, algorithms, graphs, outputDir)

	#Execute the tasks
	taskResults = map(executeTask, tasks)

	#Output files...
	writeResults(taskResults)

if __name__ == "__main__":
	sys.exit(main())
