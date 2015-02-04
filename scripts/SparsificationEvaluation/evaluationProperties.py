from scipy.spatial import distance
from networkit import *
from scipy import stats
import numpy as np
import math

# For comparison of sparsification methods, we want to calculate distance measures between (sparsified) graphs
# and properties of these graphs. Calculation methods are defined here.

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

			#Modularity
			modularity = community.Modularity().getQuality(communitiesSparsifiedGraph, sparsifiedGraph)
		else:
			randMeasure = 0.0
			nmi = 0.0
			numCommunities = 0
			minCommunitySize = 0
			maxCommunitySize = 0
			avgCommunitySize = 0
			modularity = 0.0
		return {'randMeasure':randMeasure, 'nmi':nmi,
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
		#Precalculations used below
		dd_sparsified = properties.degreeDistribution(sparsifiedGraph)
		ds_original = properties.degreeSequence(graph)
		ds_sparsified = properties.degreeSequence(sparsifiedGraph)

		#Coefficient of the sparsified graph
		fit = properties.powerlaw.Fit(dd_sparsified)
		degreeDistCoefficient = fit.alpha
		powerLawFit = properties.degreePowerLaw(sparsifiedGraph, dd_sparsified)[1]

		#Spearmans rho
		spearman_rho, spearman_p = stats.spearmanr(ds_original, ds_sparsified)

		#Relative rank error
		ranking_original = [(n, ds_original[n]) for n in graph.nodes()]
		ranking_sparsified = [(n, ds_sparsified[n]) for n in sparsifiedGraph.nodes()]
		relRankError = centrality.relativeRankError(ranking_original, ranking_sparsified)

		#Normalized absolute difference
		normalizedAbsDiff = sum([abs(ds_original[n] - ds_sparsified[n]) for n in graph.nodes()]) / graph.numberOfNodes()

		#KS D-Statistics
		ks_d, ks_p = stats.ks_2samp(ds_original, ds_sparsified)

		return {'dd_distCoefficient':degreeDistCoefficient, 'dd_powerLawFit':powerLawFit,
			'dd_spearman_rho':spearman_rho, 'dd_spearman_p':spearman_p,
			'dd_relRankError':relRankError, 'dd_normalizedAbsDiff':normalizedAbsDiff,
			'dd_ks_d':ks_d, 'dd_ks_p':ks_p}

	def getTypes(self):
		return {'dd_distCoefficient':'real', 'dd_powerLawFit':'real',
			'dd_spearman_rho':'real', 'dd_spearman_p':'real',
			'dd_relRankError':'real', 'dd_normalizedAbsDiff':'real',
			'dd_ks_d':'real', 'dd_ks_p':'real'}

#Clustering Coefficients
class P_ClusteringCoefficients:
	def getName(self):
		return "Clustering Coefficients"

	def getCCSequences(self, inputGraph):
		"""
		Returns a 2-tuple containing the list of local clustering coefficients and the list
		of average local clustering coefficients per degree.
		"""
		localCCs = properties.ClusteringCoefficient.exactLocal(inputGraph)
		maxDegree = max([inputGraph.degree(n) for n in inputGraph.nodes()])

		#Calculate average ccs per degree.
		ccDictPerDegree = {}
		for n in inputGraph.nodes():
			deg = inputGraph.degree(n)
			if not deg in ccDictPerDegree:
				ccDictPerDegree[deg] = [localCCs[n]]
			else:
				ccDictPerDegree[deg].append(localCCs[n])

		ccPerDegree = [0] * (maxDegree + 1)
		for deg in ccDictPerDegree:
			ccPerDegree[deg] = np.average(ccDictPerDegree[deg])

		return localCCs, ccPerDegree

	def getValues(self, graph, sparsifiedGraph):
		#Precalculations
		localCC_original, ccPerDegree_original = self.getCCSequences(graph)
		localCC_sparsified, ccPerDegree_sparsified = self.getCCSequences(sparsifiedGraph)

		#single-value coefficients
		ccAvgLocal = properties.ClusteringCoefficient().avgLocal(sparsifiedGraph)
		if sparsifiedGraph.numberOfNodes() < 300:
			ccGlobal = properties.ClusteringCoefficient().exactGlobal(sparsifiedGraph)
		else:
			ccGlobal = properties.ClusteringCoefficient().approxGlobal(sparsifiedGraph, min(sparsifiedGraph.numberOfNodes(), 10000))

		# ------------------------ local clustering coefficients per node ------------------------------
		#spearmans rho
		perNode_spearman_rho, perNode_spearman_p = stats.spearmanr(localCC_original, localCC_sparsified)

		#Relative rank error
		ranking_original = [(n, localCC_original[n]) for n in graph.nodes()]
		ranking_sparsified = [(n, localCC_sparsified[n]) for n in sparsifiedGraph.nodes()]
		perNode_relRankError = centrality.relativeRankError(ranking_original, ranking_sparsified)

		#Normalized absolute difference
		perNode_normalizedAbsDiff = sum([abs(localCC_original[n] - localCC_sparsified[n]) for n in graph.nodes()]) / graph.numberOfNodes()

		#KS D-Statistics
		perNode_ks_d, perNode_ks_p = stats.ks_2samp(localCC_original, localCC_sparsified)

		# ------------------------ local clustering coefficients per degree ------------------------------
		#KS D-Statistics
		perDegree_ks_d, perDegree_ks_p = stats.ks_2samp(ccPerDegree_original, ccPerDegree_sparsified)

		return {'ccAvgLocal':ccAvgLocal, 'ccGlobal':ccGlobal, 'cc_spearman_rho':perNode_spearman_rho,
			'cc_spearman_p':perNode_spearman_p, 'cc_relRankError':perNode_relRankError,
			'cc_normalizedAbsDiff':perNode_normalizedAbsDiff, 'cc_ks_d':perNode_ks_d,
			'cc_ks_p':perNode_ks_p, 'cc_perDegree_ks_d':perDegree_ks_d, 'cc_perDegree_ks_p':perDegree_ks_p }

	def getTypes(self):
		return {'ccAvgLocal':'real', 'ccGlobal':'real', 'cc_spearman_rho':'real', 'cc_spearman_p':'real',
		'cc_relRankError':'real', 'cc_normalizedAbsDiff':'real', 'cc_ks_d':'real', 'cc_ks_p':'real',
		'cc_perDegree_ks_d':'real', 'cc_perDegree_ks_p':'real' }

#Weakly Connected Components
class P_ConnectedComponents:
	def getName(self):
		return "Connected Components"

	def getComponentSizes(self, wcc):
		return list(map(lambda tuple_ID_Size: tuple_ID_Size[1], list(wcc.getComponentSizes().items())))

	def getValues(self, graph, sparsifiedGraph):
		wccs_original = properties.ConnectedComponents(graph)
		wccs_original.run()
		componentSizes_original = self.getComponentSizes(wccs_original)

		wccs_sparsified = properties.ConnectedComponents(sparsifiedGraph)
		wccs_sparsified.run()
		componentSizes_sparsified = self.getComponentSizes(wccs_sparsified)

		#Number of weakly connected components
		wccCount = wccs_sparsified.numberOfComponents()

		#KS D-Statistics
		ks_wccSizes, p_wccSizes = stats.ks_2samp(componentSizes_original, componentSizes_sparsified)

		#NMI (difference between partitions)
		nmi = community.NMIDistance().getDissimilarity(graph, wccs_original.getPartition(), wccs_sparsified.getPartition())

		return {'wcc_sizes_ks':ks_wccSizes, 'wcc_sizes_p':p_wccSizes, 'wcc_count':wccCount, 'wcc_nmi':nmi}

	def getTypes(self):
		return {'wcc_sizes_ks':'real', 'wcc_sizes_p':'real', 'wcc_count':'real', 'wcc_nmi':'real'}

#PageRank
class P_PageRank:
	def getName(self):
		return "PageRank"

	def getRanking(self, inputGraph):
		bc = centrality.PageRank(inputGraph)
		bc.run()
		return bc.ranking()

	def getValues(self, graph, sparsifiedGraph):
		ranking_original = self.getRanking(graph)
		ranking_sparsified = self.getRanking(sparsifiedGraph)

		scores_original = [r[1] for r in ranking_original]
		scores_sparsified = [r[1] for r in ranking_sparsified]
		spearman_rho, spearman_p = stats.spearmanr(scores_original, scores_sparsified)

		#Relative rank error
		relRankError = centrality.relativeRankError(ranking_original, ranking_sparsified)

		return {'pagerank_spearman_rho':spearman_rho, 'pagerank_spearman_p':spearman_p,
			'pagerank_relRankError':relRankError}

	def getTypes(self):
		return {'pagerank_spearman_rho':'real', 'pagerank_spearman_p':'real',
			'pagerank_relRankError':'real'}


# The following is not being used anymore
#Centrality
class P_Centrality_DEPRECATED:
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
		bc = centrality.ApproxBetweenness(graph, epsilon=0.05, delta=0.05, diameterSamples=0)
		#bc = centrality.Betweenness(graph)
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
class P_Components_DEPRECATED:
	def getName(self):
		return "Connected Components"

	def getValues(self, graph, sparsifiedGraph):
		nComponents, componentSizes = properties.components(sparsifiedGraph)

		return {'largestComponentSize':max(componentSizes.values()), 'numComponents':nComponents}

	def getTypes(self):
		return {'largestComponentSize':'integer', 'numComponents':'integer'}
