from scipy.spatial import distance
from scipy.stats import kendalltau
from networkit import *
import math

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
		return {'graph':'text', 'algorithm':'text', 'parameter':'text', 'evalExpr':'text', 'rt_attribute':'real', 'rt_backbone':'real'}

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

			#Number of communities
			communitySizes = communitiesBackbone.subsetSizes()
			numCommunities = communitiesBackbone.numberOfSubsets()
			minCommunitySize = min(communitySizes)
			maxCommunitySize = max(communitySizes)
			avgCommunitySize = sum(communitySizes) / len(communitySizes)

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
			numCommunities = 0
			minCommunitySize = 0
			maxCommunitySize = 0
			avgCommunitySize = 0
		return {'randMeasure':randMeasure, 'nmi':nmi, 'ccAvgLocal':ccAvgLocal, 'ccGlobal':ccGlobal,
			'numCommunities':numCommunities, 'minCommunitySize':minCommunitySize,
			'maxCommunitySize':maxCommunitySize, 'avgCommunitySize':avgCommunitySize}

	def getTypes(self):
		return {'randMeasure':'real', 'nmi':'real', 'ccAvgLocal':'real', 'ccGlobal':'real',
			'numCommunities':'integer', 'minCommunitySize':'integer',
			'maxCommunitySize':'integer', 'avgCommunitySize':'integer'
		}

#Diameter
class P_Diameter:
	def getName(self):
		return "Diameter"

	#def getValues(self, graph, backbone):
	#	diameter = properties.Diameter.estimatedDiameterRange(backbone, error=0.1)
	#	return {'diameterLow':diameter[0], 'diameterHigh':diameter[1]}

	#def getTypes(self):
	#	return {'diameterLow':'integer', 'diameterHigh':'integer'}

	def getValues(self, graph, backbone):
		diameter = properties.Diameter.estimatedDiameterRange(backbone, 0)[0]
		return {'diameter':diameter}

	def getTypes(self):
		return {'diameter':'integer'}

#Degree distribution
class P_DegreeDistribution:
	def getName(self):
		return "Degree Distribution"

	def getValues(self, graph, backbone):
		dd = properties.degreeDistribution(backbone)
		fit = properties.powerlaw.Fit(dd)
		degreeDistCoefficient = fit.alpha
		powerLawFit = properties.powerLawFit(backbone, dd)[1]

		return {'degreeDistCoefficient':degreeDistCoefficient, 'powerLawFit':powerLawFit}

	def getTypes(self):
		return {'degreeDistCoefficient':'real', 'powerLawFit':'real'}

#Centrality
class P_Centrality:
	def getName(self):
		return "Centrality"

	def getHubsFromRanking(self, ranking, count):
		ranking.sort(key=lambda x: (x[1] if not math.isnan(x[1]) else 0), reverse=True) #Sort by centrality score
		ranking = ranking[:count]
		return list(map(lambda x: x[0], ranking))

	def getBetweennessHubs(self, graph, count):
		#Empty graphs result in crash of approxbetweenness. #TODO incestigate
		if graph.numberOfNodes() == 0:
			return [-1] * count		

		print("ApproxBetweenness...")
		bc = centrality.ApproxBetweenness(graph, epsilon=0.75, delta=0.75, diameterSamples=0)
		bc.run()
		return self.getHubsFromRanking(bc.ranking(), count)

	def getPageRankHubs(self, graph, count):
		print("PageRank...")
		bc = centrality.PageRank(graph, damp=0.95)
		bc.run()
		return self.getHubsFromRanking(bc.ranking(), count)

	def getJaccard(self, list1, list2):
		return len(set(list1) & set(list2)) / len(set(list1) | set(list2))

	def getValues(self, graph, backbone):
		#lcGraph = workflows.extractLargestComponent(graph)
		#lcBackbone = workflows.extractLargestComponent(backbone)

		#PageRank
		hubCountPageRank = math.ceil(graph.numberOfNodes() * 0.01)
		prHubsG = self.getPageRankHubs(graph, hubCountPageRank)
		prHubsB = self.getPageRankHubs(backbone, hubCountPageRank)
		centralityPageRank = self.getJaccard(prHubsG, prHubsB)

		#Betweenness
		hubCountBetweenness = math.ceil(graph.numberOfNodes() * 0.001)
		bHubsG = self.getBetweennessHubs(graph, hubCountBetweenness)
		bHubsB = self.getBetweennessHubs(backbone, hubCountBetweenness)
		centralityBetweenness = self.getJaccard(bHubsG, bHubsB)

		return {'centralityPageRank':centralityPageRank, 'centralityBetweenness':centralityBetweenness}

	def getTypes(self):
		return {'centralityPageRank':'real', 'centralityBetweenness':'real'}
