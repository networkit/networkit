from NetworKit import *
from evaluation import *
from bbalgorithms import *

ABS_ZERO = 1e-5
BINARY_SEARCH_STEPS = 20

#Parameterizes the given algorithm such that the targetEdgeRatio is met approximately.
def parameterize(graph, algorithm, targetEdgeRatio):
	if algorithm.parameterizationType() == "Float":
		bestParameter = parameterizeUsingDouble(graph, algorithm, targetEdgeRatio, 0.0, 1.0, BINARY_SEARCH_STEPS)
	elif algorithm.parameterizationType() == "Int":
		bestParameter = parameterizeUsingInteger(graph, algorithm, targetEdgeRatio, 0, 10)
	elif algorithm.parameterizationType() == "Trivial":
		bestParameter = parameterizeUsingDoubleTrivial(graph, algorithm, targetEdgeRatio)
	else:
		bestParameter = None
	return bestParameter

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

		print ("Est:", str(estimation), "er: ", currentEdgeRatio, ", ter: ", targetEdgeRatio)
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

	return bestParameter

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

		print ("Current er: ", currentEdgeRatio, ", target er: ", targetEdgeRatio)

		distance = abs(currentEdgeRatio - targetEdgeRatio)
		if distance < minDistance and abs(currentEdgeRatio) > ABS_ZERO:
			minDistance = distance
			bestParameter = i
			bestRatio = currentEdgeRatio

	print("Best fit for parameter: ", bestParameter, " er: ", bestRatio, ", ter: ", targetEdgeRatio)
	return bestParameter

#Directly use targetEdgeRatio as input parameter for the algorithm.
def parameterizeUsingDoubleTrivial(G, algorithm, targetEdgeRatio):
	print("Parameterizing ", algorithm.getName(), " ... (trivial), edge ratio ", str(targetEdgeRatio))
	return targetEdgeRatio
