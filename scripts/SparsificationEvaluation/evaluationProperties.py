from scipy.spatial import distance
from scipy.stats import kendalltau
from networkit import *
from scipy import stats # for KS statistic
import numpy as np
import math

#This file holds the definitions of all backbone properties, including
#how to calculate the property's characteristic values.

#Base class
class P_Example:
	def getName(self):
		return "Example property calculator"

	#Returns a dictionary containing key/value pairs that are calculated from the given graph and sparsified graph.
	def getValues(self, graph, sparsifiedGraph):
		return {'key':'value'}

	#Returns a dictionary containing the types (integer or real) of the values
	def getTypes(self):
		return {'key':'real'}

#Fake class, ensuring some fields are created in the sqlite database
class P_General:
	def getName(self):
		return "General information"

	def getValues(self, graph, sparsifiedGraph):
		return {}

	def getTypes(self):
		return {'graph':'text', 'algorithm':'text', 'parameter':'text', 'evalExpr':'text', 'rt_attribute':'real', 'rt_backbone':'real', 'targetEdgeRatio':'real'}

#Node and edge ratios
class P_Ratios:
	def getName(self):
		return "Node and edge ratios"

	def getValues(self, graph, sparsifiedGraph):
		numNodes = sparsifiedGraph.numberOfNodes()
		numEdges = sparsifiedGraph.numberOfEdges()
		nodeRatio = (numNodes / graph.numberOfNodes())
		edgeRatio = (numEdges / graph.numberOfEdges())
		return {'numNodes':numNodes, 'numEdges':numEdges, 'nodeRatio':nodeRatio, 'edgeRatio':edgeRatio}

	def getTypes(self):
		return {'numNodes':'integer', 'numEdges':'integer', 'nodeRatio':'real', 'edgeRatio':'real'}

#Community structure
class P_Community:
	def getName(self):
		return "Community structure"

	def getValues(self, graph, sparsifiedGraph):
		if sparsifiedGraph.numberOfEdges() > 0:
			cAlgo = community.PLM(graph, refine=False, par='none')
			communitiesGraph = community.detectCommunities(graph, algo=cAlgo)
			communitiesSparsifiedGraph = community.detectCommunities(sparsifiedGraph, algo=cAlgo)

			#Number of communities
			communitySizes = communitiesSparsifiedGraph.subsetSizes()
			numCommunities = communitiesSparsifiedGraph.numberOfSubsets()
			minCommunitySize = min(communitySizes)
			maxCommunitySize = max(communitySizes)
			avgCommunitySize = sum(communitySizes) / len(communitySizes)

			#Graph structural rand measure
			_randMeasure = community.GraphStructuralRandMeasure()
			randMeasure = _randMeasure.getDissimilarity(graph, communitiesGraph, communitiesSparsifiedGraph)

			#Normalized Mutual information
			_nmi = community.NMIDistance()
			nmi = _nmi.getDissimilarity(graph, communitiesGraph, communitiesSparsifiedGraph)

			#Clustering coefficients
			_cc = properties.ClusteringCoefficient()
			ccAvgLocal = _cc.avgLocal(sparsifiedGraph)
			if sparsifiedGraph.numberOfNodes() < 300:
				ccGlobal = _cc.exactGlobal(sparsifiedGraph)
			else:
				ccGlobal = _cc.approxGlobal(sparsifiedGraph, min(sparsifiedGraph.numberOfNodes(), 10000))

			#Modularity
			modularity = community.Modularity().getQuality(communitiesSparsifiedGraph, sparsifiedGraph)
		else:
			randMeasure = 0.0
			nmi = 0.0
			ccAvgLocal = 0.0
			ccGlobal = 0.0
			numCommunities = 0
			minCommunitySize = 0
			maxCommunitySize = 0
			avgCommunitySize = 0
			modularity = 0.0
		return {'randMeasure':randMeasure, 'nmi':nmi, 'ccAvgLocal':ccAvgLocal, 'ccGlobal':ccGlobal,
			'numCommunities':numCommunities, 'minCommunitySize':minCommunitySize,
			'maxCommunitySize':maxCommunitySize, 'avgCommunitySize':avgCommunitySize,
			'modularity':modularity}

	def getTypes(self):
		return {'randMeasure':'real', 'nmi':'real', 'ccAvgLocal':'real', 'ccGlobal':'real',
			'numCommunities':'integer', 'minCommunitySize':'integer',
			'maxCommunitySize':'integer', 'avgCommunitySize':'integer',
			'modularity':'real'
		}

#Diameter
class P_Diameter:
	def getName(self):
		return "Diameter"

	#def getValues(self, graph, sparsifiedGraph):
	#	diameter = properties.Diameter.estimatedDiameterRange(sparsifiedGraph, error=0.1)
	#	return {'diameterLow':diameter[0], 'diameterHigh':diameter[1]}

	#def getTypes(self):
	#	return {'diameterLow':'integer', 'diameterHigh':'integer'}

	def getValues(self, graph, sparsifiedGraph):
		diameter = properties.Diameter.estimatedDiameterRange(sparsifiedGraph, 0)[0]
		return {'diameter':diameter}

	def getTypes(self):
		return {'diameter':'integer'}

#Degree distribution
class P_DegreeDistribution:
	def getName(self):
		return "Degree Distribution"

	def getValues(self, graph, sparsifiedGraph):
		dd = properties.degreeDistribution(sparsifiedGraph)
		fit = properties.powerlaw.Fit(dd)
		degreeDistCoefficient = fit.alpha
		powerLawFit = properties.degreePowerLaw(sparsifiedGraph, dd)[1]

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
		#bc = centrality.ApproxBetweenness(graph, epsilon=0.01, delta=0.01, diameterSamples=0)
		bc = centrality.Betweenness(graph)
		bc.run()
		return self.getHubsFromRanking(bc.ranking(), count)

	def getPageRankHubs(self, graph, count):
		print("PageRank...")
		bc = centrality.PageRank(graph, damp=0.95)
		bc.run()
		return self.getHubsFromRanking(bc.ranking(), count)

	def getJaccard(self, list1, list2):
		return len(set(list1) & set(list2)) / len(set(list1) | set(list2))

	def getValues(self, graph, sparsifiedGraph):
		#lcGraph = workflows.extractLargestComponent(graph)
		#lcSparsifiedGraph = workflows.extractLargestComponent(sparsifiedGraph)

		#PageRank
		hubCountPageRank = math.ceil(graph.numberOfNodes() * 0.01)
		prHubsG = self.getPageRankHubs(graph, hubCountPageRank)
		prHubsB = self.getPageRankHubs(sparsifiedGraph, hubCountPageRank)
		centralityPageRank = self.getJaccard(prHubsG, prHubsB)

		#Betweenness
		hubCountBetweenness = math.ceil(graph.numberOfNodes() * 0.001)
		bHubsG = self.getBetweennessHubs(graph, hubCountBetweenness)
		bHubsB = self.getBetweennessHubs(sparsifiedGraph, hubCountBetweenness)
		centralityBetweenness = self.getJaccard(bHubsG, bHubsB)

		return {'centralityPageRank':centralityPageRank, 'centralityBetweenness':centralityBetweenness}

	def getTypes(self):
		return {'centralityPageRank':'real', 'centralityBetweenness':'real'}

#Connected components
class P_Components:
	def getName(self):
		return "Connected Components"

	def getValues(self, graph, sparsifiedGraph):
		nComponents, componentSizes = properties.components(sparsifiedGraph)

		return {'largestComponentSize':max(componentSizes.values()), 'numComponents':nComponents}

	def getTypes(self):
		return {'largestComponentSize':'integer', 'numComponents':'integer'}

#Various properties based on the Kolmogorow-Smirnow-Test
class P_KolmogorowSmirnow:
	def getName(self):
		return "KolmogorowSmirnow"

	def getWCCSizes(self, inputGraph):
		wccs = properties.ConnectedComponents(inputGraph)
		wccs.run()
		componentSizes = list(map(lambda tuple_ID_Size: tuple_ID_Size[1], list(wccs.getComponentSizes().items())))
		#componentSizesDist = list(map(lambda s : len([c for c in componentSizesList if c == s]), range(0, max(componentSizesList) + 1)))
		return componentSizes

	def getValues(self, graph, sparsifiedGraph):
		#Degree Distribution
		sampleGraph = properties.degreeSequence(graph)
		sampleSparsifiedGraph = properties.degreeSequence(sparsifiedGraph)
		ks_dd, p_dd = stats.ks_2samp(sampleGraph, sampleSparsifiedGraph)

		#Distribution of clustering coefficients (TODO: per degree?)
		sampleGraph = properties.ClusteringCoefficient.exactLocal(graph)#self.getCCDistPerDegree(graph)
		sampleSparsifiedGraph = properties.ClusteringCoefficient.exactLocal(sparsifiedGraph)#self.getCCDistPerDegree(sparsifiedGraph)
		ks_cc_perDegree, p_cc_perDegree = stats.ks_2samp(sampleGraph, sampleSparsifiedGraph)

		#Distribution of sizes of weakly connected components
		sampleGraph = self.getWCCSizes(graph)
		sampleSparsifiedGraph = self.getWCCSizes(sparsifiedGraph)
		ks_wccSizes, p_wccSizes = stats.ks_2samp(sampleGraph, sampleSparsifiedGraph)

		return {'ks_dd':ks_dd, 'p_dd':p_dd, 'ks_cc_perDegree':ks_cc_perDegree, 'p_cc_perDegree':p_cc_perDegree, 'ks_wccSizes':ks_wccSizes, 'p_wccSizes':p_wccSizes}

	def getTypes(self):
		return {'ks_dd':'real', 'p_dd':'real', 'ks_cc_perDegree':'real', 'p_cc_perDegree':'real', 'ks_wccSizes':'real', 'p_wccSizes':'real'}
