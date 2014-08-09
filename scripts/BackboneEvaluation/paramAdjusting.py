from NetworKit import *

ABS_ZERO = 1e-5

def parameterizeUsingDouble(G, algorithm, targetEdgePercentage, lowerBound, upperBound, maxSteps):
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
		currentEdgePercentage = backbone.numberOfEdges() / G.numberOfEdges()

		print ("Est:", str(estimation), "ep: ", currentEdgePercentage, ", tep: ", targetEdgePercentage)
		distance = abs(currentEdgePercentage - targetEdgePercentage)

		if distance < minDistance and abs(currentEdgePercentage) > ABS_ZERO:
			minDistance = distance
			bestParameter = estimation

		#"Exact" hit?
		if abs(currentEdgePercentage - targetEdgePercentage) < 1e-7:
			break;

		increase = (currentEdgePercentage < targetEdgePercentage)
		if not algorithm.increasing():
			increase = not increase

		if increase:
			lowerBound = estimation
		else:
			upperBound = estimation

	return algorithm.getAlgorithmExpr(bestParameter)


def parameterizeUsingInteger(G, algorithm, targetEdgePercentage, lowerBound, upperBound):
	print("Parameterizing ", algorithm.getName(), " ...")
	bestParameter = lowerBound
	bestPercentage = 0.0
	minDistance = 100.0

	#Precalculation
	print ("Calculating attribute...")
	preCalcAttr = algorithm.getPrecalcAttribute(G)

	print("Fitting parameter...")
	for i in range(lowerBound, upperBound + 1):
		backbone = algorithm.getPrecalcBackbone(G, preCalcAttr, i)
		currentEdgePercentage = backbone.numberOfEdges() / G.numberOfEdges()

		print ("Current ep: ", currentEdgePercentage, ", target ep: ", targetEdgePercentage)

		distance = abs(currentEdgePercentage - targetEdgePercentage)
		if distance < minDistance and abs(currentEdgePercentage) > ABS_ZERO:
			minDistance = distance
			bestParameter = i
			bestPercentage = currentEdgePercentage

	print("Best fit for parameter: ", bestParameter, " ep: ", bestPercentage, ", tep: ", targetEdgePercentage)
	return algorithm.getAlgorithmExpr(bestParameter)
