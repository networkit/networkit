from NetworKit import *
from scipy.spatial import distance
import time

# -----------------------------------------------------------------------
# The purpose of the following script is to automatically apply a set
# of backbone algorithms to a set of input graphs and to determine
# certain graph properties for evaluation.
# -----------------------------------------------------------------------


#A pair of a graph description and a set of backbone algorithms
class BackboneTask:
	def __init__(self, graph, algorithms, outputDir=None):
		self.graph = graph
		self.algorithms = algorithms
		self.outputDir = outputDir

class BackboneTaskResult:
	def __init__(self, task):
		self.loadingTime = 0
		self.backboneProperties = []
		self.task = task

#Information about a graph; used as input parameter
class GraphDescription:
	def __init__(self, path, format, name):
		self.path = path
		self.format = format
		self.name = name

class BackboneAlgorithm:
	def __init__(self, algorithmString, name):
		self.algorithmString = algorithmString
		self.name = name

#Contains some result values that are characteristic for a backbone calculation result
class GraphProperties:
	def __init__(self, algorithmName):
		self.numNodes = 0
		self.keptNodesPercent = 0
		self.numEdges = 0
		self.keptEdgesPercent = 0
		self.graphStructuralRandMeasure = 0
		self.nmi = 0
		self.clusteringCoefficient = 0
		self.degreeDistCoefficient = 0
		self.cpvDistanceFromOriginal = 0
		self.cpvDistanceFromOriginalNormalized = 0
		self.runningTime = 0.0
		self.diameter = None
		self.algorithmName = algorithmName

# Fake backbone algorithm that returns the input graph itself
class OriginalGraph:
	def calculate(self, graph):
		return graph

def debugInfo(text):
	print("{" + text + "...}")

def getCentralityPositionVector(graph):
	bc = centrality.ApproxBetweenness2(graph, min(100, graph.numberOfNodes()))
	bc.run()
	ranking = map(lambda x: x[0], bc.ranking())
	centralityPositionVector = [0] * graph.numberOfNodes()
	rank = 0
	for node in ranking:
		centralityPositionVector[node] = rank
		rank += 1
	return centralityPositionVector

# Applies the given backbone algorithm to the given graph and calculates various graph properties
def applyBackboneAlgorithm(graph, algorithm):
	#Backbone calculation
	debugInfo("Calculating backbone")
	start = time.clock()
	backbone = eval(algorithm.algorithmString).calculate(graph)
	end = time.clock()

	#Result backbone graph properties
	bprops = GraphProperties(algorithm.name)
	bprops.runningTime = end - start

	#Basic graph properties
	bprops.numNodes = backbone.numberOfNodes()
	bprops.numEdges = backbone.numberOfEdges()
	bprops.keptNodesPercent = (bprops.numNodes / graph.numberOfNodes()) * 100.0
	bprops.keptEdgesPercent = (bprops.numEdges / graph.numberOfEdges()) * 100.0
	debugInfo("Backbone edge percentage: " + str(bprops.keptEdgesPercent))

	debugInfo("Detecting communities")
	if backbone.numberOfEdges() > 0:
		communitiesGraph = community.detectCommunities(graph)
		communitiesBackbone = community.detectCommunities(backbone)

		#Graph structural rand measure
		randMeasure = community.GraphStructuralRandMeasure()
		bprops.graphStructuralRandMeasure = randMeasure.getDissimilarity(graph, communitiesGraph, communitiesBackbone)

		#Normalized Mutual information
		nmi = community.NMIDistance()
		bprops.nmi = nmi.getDissimilarity(graph, communitiesGraph, communitiesBackbone)

		#Clustering coefficient
		debugInfo("Calculating clustering coefficient")
		#bprops.clusteringCoefficient = properties.clustering(backbone)
		cc = properties.ClusteringCoefficient()
		bprops.clusteringCoefficient = cc.avgLocal(backbone)

	#Diameter
	debugInfo("Calculating diameter")
	bprops.diameter = str(properties.Diameter.estimatedDiameterRange(backbone, error=0.1))

	#Degree distribution coefficient
	debugInfo("Calculating degree distribution coefficient")
	bprops.degreeDistCoefficient = properties.powerLawExponent(backbone)

	#Centrality position vector distance
	debugInfo("Calculating centrality position vector distance")
	cpvOriginal = getCentralityPositionVector(graph)
	cpvBackbone = getCentralityPositionVector(backbone)
	bprops.cpvDistanceFromOriginal = distance.euclidean(cpvOriginal, cpvBackbone)
	bprops.cpvDistanceFromOriginalNormalized = bprops.cpvDistanceFromOriginal / graph.numberOfNodes()

	return backbone, bprops

def writeBackboneFile(backbone, task, algorithm):
	fileName = task.outputDir + "/" + task.graph.name + "_" + algorithm.name + ".txt"
	if task.outputDir is not None:
		writeGraph(backbone, fileName, Format.METIS)


# Executes the given BackboneTask
def executeTask(task):
	#Load graph into memory
	start = time.clock()
	graph = readGraph(task.graph.path, task.graph.format)
	end = time.clock()

	taskResult = BackboneTaskResult(task)
	taskResult.loadingTime = end - start

	#Apply the algorithms!
	for algorithm in task.algorithms:
		print("[Applying algorithm '", algorithm.name, "' to graph '", task.graph.name , "']")
		backbone, bprops = applyBackboneAlgorithm(graph, algorithm)
		taskResult.backboneProperties.append(bprops)

		#Write output
		writeBackboneFile(backbone, task, algorithm)

	return taskResult
