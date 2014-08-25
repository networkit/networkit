from scipy.spatial import distance
from scipy.stats import kendalltau
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
        if backbone.numberOfEdges() > 0:
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
        else:
            randMeasure = 0.0
            nmi = 0.0
            ccAvgLocal = 0.0
            ccGlobal = 0.0
        return {'randMeasure':randMeasure, 'nmi':nmi, 'ccAvgLocal':ccAvgLocal, 'ccGlobal':ccGlobal}

    def getTypes(self):
        return {'randMeasure':'real', 'nmi':'real', 'ccAvgLocal':'real', 'ccGlobal':'real'}


#Diameter
class P_Diameter:
    def getName(self):
        return "Diameter"

    #def getValues(self, graph, backbone):
    #    diameter = properties.Diameter.estimatedDiameterRange(backbone, error=0.1)
    #    return {'diameterLow':diameter[0], 'diameterHigh':diameter[1]}

    #def getTypes(self):
    #    return {'diameterLow':'integer', 'diameterHigh':'integer'}

    def getValues(self, graph, backbone):
        if graph.isWeighted():
            diameter = properties.Diameter.exactDiameter(workflows.extractLargestComponent(backbone))
        else:
            #This is actually not neccessary but a workaround for sometimes failing diameter calculation. TODO: investigate
            diameter = properties.Diameter.exactDiameter(backbone)
        return {'diameter':diameter}

    def getTypes(self):
        return {'diameter':'integer'}

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

    #Returns a list of node ids. The nodes are sorted in descending order by betweenness centrality.
    def getCentralityRanking(self, graph):
        bc = centrality.ApproxBetweenness2(graph, min(100, graph.numberOfNodes()))
        bc.run()
        ranking = bc.ranking()
        ranking.sort(key=lambda x: x[1]) #Sort by centrality score
        return list(map(lambda x: x[0], ranking))

    def kendallTauQM(self, c1, c2):
        # Transform from [-1, 1] to [0, 1]. 0 is worst and returned for reverse ordering.
        return (kendalltau(c1, c2)[0] + 1.0) / 2.0

    def getValues(self, graph, backbone):
        cpvOriginal = self.getCentralityPositionVector(graph)
        cpvBackbone = self.getCentralityPositionVector(backbone)
        cpvDistance = distance.euclidean(cpvOriginal, cpvBackbone)
        cpvDistanceNormalized = cpvDistance / graph.numberOfNodes()

        rankingOriginal = self.getCentralityRanking(graph)
        rankingBackbone = self.getCentralityRanking(backbone)
        bcKendallTau = self.kendallTauQM(rankingOriginal, rankingBackbone)

        return {'cpvDistance':cpvDistance, 'cpvDistanceNormalized':cpvDistanceNormalized, 'bcKendallTau':bcKendallTau}

    def getTypes(self):
        return {'cpvDistance':'real', 'cpvDistanceNormalized':'real', 'bcKendallTau':'real'}
