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
	estimation = lowerBound
	bestParameter = lowerBound
	minDistance = upperBound

	#Precalculation
	preCalcAttr = algorithm.getAttribute(G)

	for i in range(0, maxSteps):
		estimation = (lowerBound + upperBound) / 2.0
		backbone = algorithm.getBackboneFromAttribute(G, preCalcAttr, estimation)
		currentEdgeRatio = backbone.numberOfEdges() / G.numberOfEdges()

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
	bestParameter = lowerBound
	bestRatio = 0.0
	minDistance = 100.0

	#Precalculation
	preCalcAttr = algorithm.getAttribute(G)

	for i in range(lowerBound, upperBound + 1):
		backbone = algorithm.getBackboneFromAttribute(G, preCalcAttr, i)
		currentEdgeRatio = backbone.numberOfEdges() / G.numberOfEdges()

		distance = abs(currentEdgeRatio - targetEdgeRatio)
		if distance < minDistance and abs(currentEdgeRatio) > ABS_ZERO:
			minDistance = distance
			bestParameter = i
			bestRatio = currentEdgeRatio

	return bestParameter

#Directly use targetEdgeRatio as input parameter for the algorithm.
def parameterizeUsingDoubleTrivial(G, algorithm, targetEdgeRatio):
	return targetEdgeRatio
