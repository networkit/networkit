from networkit import *
from evaluation import *
from txtWriter import TxtResultWriter
from sqliteWriter import SqliteResultWriter
from consoleWriter import ConsoleWriter
import string
import parameters
import multiprocessing

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

def main():
	#Get input parameters
	edgeRatios = parameters.getEdgeRatios()
	algorithms = parameters.getAlgorithms()
	graphs = parameters.getGraphs()
	properties = parameters.getProperties()

	#Generate tasks from input parameters
	lock = multiprocessing.Lock()
	writers = [SqliteResultWriter("./output/backbones.db")]
	jobs = [multiprocessing.Process(target=executeTask, name=g, args=(Task([g], algorithms, properties, edgeRatios), writers, lock)) for g in graphs]

	#Execute the tasks
	for job in jobs:
		job.start()

	#Wait for the tasks to finish
	for job in jobs:
		job.join()

	print("[FINISHED ALL JOBS]")

if __name__ == "__main__":
	sys.exit(main())
