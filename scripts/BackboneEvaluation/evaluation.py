from NetworKit import *
from scipy.spatial import distance
import time

# -----------------------------------------------------------------------
# The purpose of the following script is to automatically apply a set
# of backbone algorithms to a set of input graphs and to determine
# certain graph properties for evaluation.
# -----------------------------------------------------------------------

#Sets of graphs, algorithms, properties, and target edge ratios.
class Task:
	def __init__(self, graphs, algorithms, properties, edgeRatios):
		self.graphs = graphs
		self.algorithms = algorithms
		self.properties = properties
		self.edgeRatios = edgeRatios

#Holds the results of a Task
class TaskResult:
	def __init__(self, task):
		self.task = task
		self.data = [] #Set of dictionaries

#Information about a graph; used as input parameter
class GraphDescription:
	def __init__(self, path, format, name):
		self.path = path
		self.format = format
		self.name = name

# Fake backbone algorithm that returns the input graph itself
class OriginalGraph:
	def calculate(self, graph):
		return graph


# Calculates all backbone properties for all graphs for all algorithms.
def executeTask(task):
	taskResult = TaskResult(task)
	data = [] 				#Set of dictionaries containing key/value pairs
	columns = [] 		#Set containing the keys that appear in data.

	for igraph in task.graphs:
		graph = readGraph(igraph.path, igraph.format)

		for ialgorithm in task.algorithms:
			#Calculate the attribute that is characteristic for that algorithm.
			attribute = ialgorithm.getPrecalcAttribute(graph)

			if not G.isWeighted() and inputAlgorithm.requiresWeight():
				print("Skipping ", igraph.name, " for ", ialgorithm.getName(), " (requires weighted graph)")
				continue

			for iedgeRatio in task.edgeRatios:
				#Parameterize the algorithm in such a way that we meet the expected edge ratio
				algorithmParameter = parameterization.parameterize(graph, ialgorithm, iedgeRatio)


	#Apply the algorithms!
	for algorithm in task.algorithms:
		print("[Applying algorithm '", algorithm.name, "' to graph '", task.graph.name , "']")
		backbone, bprops = applyBackboneAlgorithm(graph, algorithm)
		taskResult.backboneProperties.append(bprops)

		#Write output
		writeBackboneFile(backbone, task, algorithm)

	return taskResult
