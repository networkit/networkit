from NetworKit import *
from evaluation import *
from bbalgorithms import *

ABS_ZERO = 1e-5

# Returns backbone task that contains all combinations of desired algorithms and target edge ratios
# for the given graph.
def getBackboneTask(graph, inputAlgorithms, targetEdgeRatios, outputDir):
	#TODO: actually use those interval limits given in the input strings...

	#Load graph
	G = readGraph(graph.path, graph.format)
	G.indexEdges()

	#Generate a set of algorithms to apply to that graph
	algorithms = []
	for inputAlgorithm in inputAlgorithms:
		#Skip unweighted graphs for algorithms that require a weight
		if not G.isWeighted() and inputAlgorithm.requiresWeight():
			print("Skipping ", graph.name, " for ", inputAlgorithm.getName(), " (requires weighted graph)")
			continue

		for targetEp in targetEdgeRatios:
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
	return BackboneTask(graph, algorithms, outputDir)

#Parameterize the given algorithm using binary search.
def parameterizeUsingDouble(G, algorithm, targetEdgeRatio, lowerBound, upperBound, maxSteps):
	print("Parameterizing ", algorithm.getName(), " ...")
	estimation = lowerBound
	bestParameter = lowerBound
	minDistance = upperBound

	#Precalculation
	print ("Calculating attribute...")
	preCalcAttr = algorithm.getPrecalcAttribute(G)

	print("Fitting parameter...")
	for i in range(0, maxSteps):
		estimation = (lowerBound + upperBound) / 2.0
		backbone = algorithm.getPrecalcBackbone(G, preCalcAttr, estimation)
		currentEdgeRatio = backbone.numberOfEdges() / G.numberOfEdges()

		print ("Est:", str(estimation), "ep: ", currentEdgeRatio, ", tep: ", targetEdgeRatio)
		distance = abs(currentEdgeRatio - targetEdgeRatio)

		if distance < minDistance and abs(currentEdgeRatio) > ABS_ZERO:
			minDistance = distance
			bestParameter = estimation

		#"Exact" hit?
		if abs(currentEdgeRatio - targetEdgeRatio) < 1e-7:
			break;

		increase = (currentEdgeRatio < targetEdgeRatio)
		if not algorithm.increasing():
			increase = not increase

		if increase:
			lowerBound = estimation
		else:
			upperBound = estimation

	return algorithm.getAlgorithmExpr(bestParameter)

#Parameterize the given algorithm using brute force search.
def parameterizeUsingInteger(G, algorithm, targetEdgeRatio, lowerBound, upperBound):
	print("Parameterizing ", algorithm.getName(), " ...")
	bestParameter = lowerBound
	bestRatio = 0.0
	minDistance = 100.0

	#Precalculation
	print ("Calculating attribute...")
	preCalcAttr = algorithm.getPrecalcAttribute(G)

	print("Fitting parameter...")
	for i in range(lowerBound, upperBound + 1):
		backbone = algorithm.getPrecalcBackbone(G, preCalcAttr, i)
		currentEdgeRatio = backbone.numberOfEdges() / G.numberOfEdges()

		print ("Current ep: ", currentEdgeRatio, ", target ep: ", targetEdgeRatio)

		distance = abs(currentEdgeRatio - targetEdgeRatio)
		if distance < minDistance and abs(currentEdgeRatio) > ABS_ZERO:
			minDistance = distance
			bestParameter = i
			bestRatio = currentEdgeRatio

	print("Best fit for parameter: ", bestParameter, " ep: ", bestRatio, ", tep: ", targetEdgeRatio)
	return algorithm.getAlgorithmExpr(bestParameter)
