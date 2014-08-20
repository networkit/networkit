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

def writeResults(taskResults):
	#TxtResultWriter("./output"),
	#SqliteResultWriter("./output/backbones.db")
	writers = [ConsoleWriter()]
	for taskResult in taskResults:
		for writer in writers:
			writer.receiveResult(taskResult)

def main():
	#Get input parameters
	edgeRatios = parameters.getEdgeRatios()
	algorithms = parameters.getAlgorithms()
	graphs = parameters.getGraphs()
	properties = parameters.getProperties()

	#Generate tasks from input parameters
	task = Task(graphs, algorithms, properties, edgeRatios)

	#Execute the tasks
	taskResult = executeTask(task)

	#Output files...
	writeResults(taskResults)

if __name__ == "__main__":
	sys.exit(main())
