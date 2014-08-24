from scipy.spatial import distance
from networkit import *

#This file holds the definitions of all backbone properties, including
#how to calculate the property's characteristic values.

#Base class
class P_Example:
    def getName(self):
        return "Example property calculator"

    #Returns a dictionary containing key/value pairs that are calculated from the given graph and backbone.
    def getValues(self, graph, backbones):
        return {'key':'value'}

    #Returns a dictionary containing the types (integer or real) of the values
    def getTypes(self):
        return {'key':'real'}

#Fake class, ensuring some fields are created in the sqlite database
class P_General:
    def getName(self):
        return "General information"

    def getValues(self, graph, backbone):
        return {}

    def getTypes(self):
        return {'graph':'text', 'algorithm':'text', 'parameter':'text', 'evalExpr':'text'}

#Node and edge ratios
class P_Ratios:
    def getName(self):
        return "Node and edge ratios"

    def getValues(self, graph, backbone):
        numNodes = backbone.numberOfNodes()
        numEdges = backbone.numberOfEdges()
        nodeRatio = (numNodes / graph.numberOfNodes())
        edgeRatio = (numEdges / graph.numberOfEdges())
        return {'numNodes':numNodes, 'numEdges':numEdges, 'nodeRatio':nodeRatio, 'edgeRatio':edgeRatio}

    def getTypes(self):
        return {'numNodes':'integer', 'numEdges':'integer', 'nodeRatio':'real', 'edgeRatio':'real'}

#Community structure
class P_Community:
    def getName(self):
        return "Community structure"

    def getValues(self, graph, backbone):
        if backbone.numberOfEdges() == 0:
            raise Exception('Empty backbones are not allowed.')

        cAlgo = community.PLM(refine=False, par='none')
        communitiesGraph = community.detectCommunities(graph, algo=cAlgo)
        communitiesBackbone = community.detectCommunities(backbone, algo=cAlgo)

        #Graph structural rand measure
        _randMeasure = community.GraphStructuralRandMeasure()
        randMeasure = _randMeasure.getDissimilarity(graph, communitiesGraph, communitiesBackbone)

        #Normalized Mutual information
        _nmi = community.NMIDistance()
        nmi = _nmi.getDissimilarity(graph, communitiesGraph, communitiesBackbone)

        #Clustering coefficients
        _cc = properties.ClusteringCoefficient()
        ccAvgLocal = _cc.avgLocal(backbone)
        if backbone.numberOfNodes() < 300:
            ccGlobal = _cc.exactGlobal(backbone)
        else:
            ccGlobal = _cc.approxGlobal(backbone, min(backbone.numberOfNodes(), 10000))

        return {'randMeasure':randMeasure, 'nmi':nmi, 'ccAvgLocal':ccAvgLocal, 'ccGlobal':ccGlobal}

    def getTypes(self):
        return {'randMeasure':'real', 'nmi':'real', 'ccAvgLocal':'real', 'ccGlobal':'real'}


#Diameter
class P_Diameter:
    def getName(self):
        return "Diameter"

    def getValues(self, graph, backbone):
        diameter = properties.Diameter.estimatedDiameterRange(backbone, error=0.1)
        return {'diameterLow':diameter[0], 'diameterHigh':diameter[1]}

    def getTypes(self):
        return {'diameterLow':'integer', 'diameterHigh':'integer'}

#Degree distribution
class P_DegreeDistribution:
    def getName(self):
        return "Degree Distribution"

    def getValues(self, graph, backbone):
        degreeDistCoefficient = properties.powerLawExponent(backbone)
        return {'degreeDistCoefficient':degreeDistCoefficient}

    def getTypes(self):
        return {'degreeDistCoefficient':'real'}

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

    def getTypes(self):
        return {'cpvDistance':'real', 'cpvDistanceNormalized':'real'}
