from NetworKit import *

#This file holds the definitions of all backbone properties, including
#how to calculate the property's characteristic values.

#Base class
class P_Example:
    def getName(self):
        return "Example property calculator"

    #Returns a dictionary containing key/value pairs that are calculated from the given graph and backbone.
    def getValues(self, graph, backbones):
        return {}

#Node and edge ratios
class P_Ratios:
    def getName(self):
        return "Node and edge ratios"

    def getValues(self, graph, backbone):
        numNodes = backbone.numberOfNodes()
        numEdges = backbone.numberOfEdges()
        nodeRatio = (bprops.numNodes / graph.numberOfNodes())
        edgeRatio = (bprops.numEdges / graph.numberOfEdges())
        return {'numNodes':numNodes, 'numEdges':numEdges, 'nodeRatio':nodeRatio, 'edgeRatio':edgeRatio}

#Community structure
class P_Community:
    def getName(self):
        return "Community structure"

    def getValues(self, graph, backbone):
        if backbone.numberOfEdges() == 0:
            raise Exception('Empty backbones are not allowed.')

        communitiesGraph = community.detectCommunities(graph)
		communitiesBackbone = community.detectCommunities(backbone)

		#Graph structural rand measure
		_randMeasure = community.GraphStructuralRandMeasure()
		randMeasure = _randMeasure.getDissimilarity(graph, communitiesGraph, communitiesBackbone)

		#Normalized Mutual information
		_nmi = community.NMIDistance()
		nmi = _nmi.getDissimilarity(graph, communitiesGraph, communitiesBackbone)

		#Clustering coefficients
		_cc = properties.ClusteringCoefficient()
		ccAvgLocal = _cc.avgLocal(backbone)
        ccGlobal = _cc.approxGlobal(backbone, min(backbone.numberOfNodes(), 10000))

        return {'randMeasure':randMeasure, 'nmi':nmi, 'ccAvgLocal':ccAvgLocal, 'ccGlobal':ccGlobal}

#Diameter
class P_Diameter:
    def getName(self):
        return "Diameter"

    def getValues(self, graph, backbone):
        diameter = properties.Diameter.estimatedDiameterRange(backbone, error=0.1)
        return {'diameterLow':diameter[0], 'diamterHigh':diameter[1]}

#Degree distribution
class P_DegreeDistribution:
    def getName(self):
        return "Degree Distribution"

    def getValues(self, graph, backbone):
        degreeDistCoefficient = properties.powerLawExponent(backbone)
        return {'degreeDistCoefficient':degreeDistCoefficient}

#Centrality
class P_Centrality:
    def getName(self):
        return "Centrality"

    def getCentralityPositionVector(self, graph):
        bc = centrality.ApproxBetweenness2(graph, min(100, graph.numberOfNodes()))
        bc.run()
        ranking = map(lambda x: x[0], bc.ranking())
        centralityPositionVector = [0] * graph.numberOfNodes()
        rank = 0
        for node in ranking:
            centralityPositionVector[node] = rank
            rank += 1
        return centralityPositionVector

    def getValues(self, graph, backbone):
        cpvOriginal = self.getCentralityPositionVector(graph)
        cpvBackbone = self.getCentralityPositionVector(backbone)
        cpvDistance = distance.euclidean(cpvOriginal, cpvBackbone)
        cpvDistanceNormalized = cpvDistance / graph.numberOfNodes()

        return {'cpvDistance':cpvDistance, 'cpvDistanceNormalized':cpvDistanceNormalized}
