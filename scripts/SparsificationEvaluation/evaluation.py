from networkit import *
import time
import traceback

# -----------------------------------------------------------------------
# The purpose of the following script is to automatically apply a set
# of sparsification algorithms to a set of input graphs and to determine
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
		self.data = [] 				#Set of dictionaries containing key/value pairs
		self.columns = []		#Set containing the keys that appear in data. #TODO

#Information about a graph; used as inpu tparamete
class GraphDescription:
	def __init__(self, path, format, name):
		self.path = path
		self.format = format
		self.name = name

# Calculates all properties for all graphs for all algorithms.
def executeTask(task):
	taskResult = TaskResult(task)
	taskResult.data = []
	taskResult.columns = []

	for igraph in task.graphs:
		print("-------------------------------------------------------------------------------------------- ")
		print("Now working on graph ", igraph.name)
		print("-------------------------------------------------------------------------------------------- ")
		graph = readGraph(igraph.path, igraph.format)
		graph.indexEdges()

		for ialgorithm in task.algorithms:
			#Calculate the attribute that is characteristic for that algorithm.
			time_attribute_start = time.time()
			attribute = ialgorithm.getAlgorithm().getAttribute(graph)
			time_attribute_elapsed = time.time() - time_attribute_start

			#Check preconditions
			if not graph.isWeighted() and ialgorithm.requiresWeight():
				print("Skipping ", igraph.name, " for ", ialgorithm.getShortName(), " (requires weighted graph)")
			else:
				for iedgeRatio in task.edgeRatios:
					#Parameterize the algorithm in such a way that we meet the expected edge ratio
					algorithmParameter = ialgorithm.getAlgorithm().getParameter(graph, iedgeRatio)
					
					time_backbone_start = time.time()
					backbone = ialgorithm.getAlgorithm().getSparsifiedGraph(graph, algorithmParameter, attribute)
					time_backbone_elapsed = time.time() - time_backbone_start

					propertiesDict = {
								'graph':igraph.name,
								'algorithm':ialgorithm.getShortName(),
								'parameter':algorithmParameter,
								'evalExpr':'unknown',
								'rt_attribute':time_attribute_elapsed,
								'rt_backbone':time_backbone_elapsed,
								'targetEdgeRatio':iedgeRatio
								}

					#Calculate all desired properties of the backbone
					print('Calculating properties: ', igraph.name, ', ', ialgorithm.getShortName(),', ', algorithmParameter)
					for iproperty in task.properties:
						try:
							d = iproperty.getValues(graph, backbone)
							propertiesDict = dict(list(propertiesDict.items()) + list(d.items()))
						except:
							print(traceback.format_exc())
							print("Unexpected error:", sys.exc_info()[0])
							raise

					taskResult.data.append(propertiesDict)

					#Non-parameterizable algorithms do not need to be applied several times.
					if ialgorithm.parameterizationType() == 'None':
						break

	return taskResult
