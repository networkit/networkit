"""
This module provides graph generators that produce synthetic networks according to various models.
"""

__author__ = "Christian Staudt"

# extension imports
from _NetworKit import Graph, BarabasiAlbertGenerator, PubWebGenerator, ErdosRenyiGenerator, ClusteredRandomGraphGenerator, DorogovtsevMendesGenerator, DynamicPubWebGenerator, DynamicPathGenerator, ChungLuGenerator, HyperbolicGenerator, DynamicHyperbolicGenerator, HavelHakimiGenerator, DynamicDorogovtsevMendesGenerator, RmatGenerator, DynamicForestFireGenerator, RegularRingLatticeGenerator, WattsStrogatzGenerator, PowerlawDegreeSequence, EdgeSwitchingMarkovChainGenerator, EdgeSwitchingMarkovChainGenerator as ConfigurationModelGenerator, LFRGenerator

from networkit import distance, coarsening, matching, nxadapter, graphio, graph

import math
import logging
import subprocess
import os

class StarGraphGenerator:
	"""
	Creates a graph with star topology, i.e. one
	central node connected to all other nodes.
	"""

	def __init__(self, n):
		self.n = n

	def generate(self):
		G = Graph(self.n)
		for u in G.nodes():
			if u is 0:
				pass
			else:
				G.addEdge(u, 0)
		return G


class MultiscaleGenerator:
	"""
	TODO:
	"""

	def __init__(self, O, limitLevels=4, withADWeights=True):
		self.O = O				# original graph

		# hierarchy
		self.Gc = []			# hierarchy of coarse graphs: self.Gc[self.maxLevel]
		self.Gf = []			# hierarchy of fine graphs
		self.matching = []
		self.aggregates = []	# hierarchy of aggregates
		self.up = []			# (in v cycle): up[i][u] maps node u of fine graph at level i to set of nodes of finer graph at level i+1
		self.down = []			# (in v cycle): down[i][u] maps node u of coarse graph at level i to node v of coarser graph at level i+1
		self.nodeWeights = []

		# state
		self.maxLevel = 0		# actual index of highest level

		# parameters
		self.withADWeights = True
		self.limitLevels = limitLevels
		self.aggregationScheme = "matching"


	def _weightsFromADScores(self, scores, G):
		""" turns (modified) AD distance scores into edge weights for matching """
		epsilon = 1e-12
		weights = [None for i in range(G.upperEdgeIdBound())]
		def setWeight(u,v,w,eid):
			weights[eid] = (1 / (scores[eid] + epsilon) * math.sqrt(G.degree(u) * G.degree(v)))
		G.forEdges(lambda u,v,w,eid: setWeight(u,v,w,eid))
		return weights

	def _printHierarchy(self):
		""" DEBUG method: print the coarsening hierarchy"""
		for i in range(self.maxLevel):
			print("{0}\t\t{1}\t\t{2}".format(i, self.Gc[i], self.Gf[i]))

	def _showCoarseSequence(self):
		import matplotlib.pyplot as plt
		import matplotlib.gridspec as gridspec
		import networkx

		k = 4
		gs = gridspec.GridSpec(k, k)
		plt.figure(figsize=(12,12))
		for i in range(self.maxLevel+1):
			plt.subplot(gs[math.floor(i / k), i % k])
			networkx.draw(nxadapter.nk2nx(self.Gc[i]), node_size=[(self.nodeWeights[i][v] + 5) for v in self.Gc[i].nodes()], node_color="gray")

	def _showFineSequence(self):
		import matplotlib.pyplot as plt
		import matplotlib.gridspec as gridspec
		import networkx

		k = 4
		gs = gridspec.GridSpec(k, k)
		plt.figure(figsize=(12,12))
		for i in range(self.maxLevel, -1, -1):
			print("level ", i)
			plt.subplot(gs[math.floor(i / k), i % k])
			nodeSizes = [(self.nodeWeights[i][v] + 5) for v in self.Gf[i].nodes()]
			networkx.draw(nxadapter.nk2nx(self.Gf[i]), node_size=nodeSizes, node_color="gray")

	def _buildCoarseSequence(self):
		logging.info("building coarse sequence")
		for i in range(self.limitLevels):
			self.maxLevel += 1
			logging.info("level: ", i)
			if i is 0:
				self.Gc.append(self.O)
				self.down.append({})
				self.up.append({})
				self.nodeWeights.append([1 for v in range(self.Gc[i].upperNodeIdBound())])
			else:
				if self.aggregationScheme == "matching":
					# perform matching
					logging.info("\t matching")
					if self.withADWeights:
						logging.info("\t calculating algebraic distance weights")
						# index edges if not already happened
						if not self.Gc[i-1].hasEdgeIds():
							self.Gc[i-1].indexEdges()
						ad = distance.AlgebraicDistance(self.Gc[i-1], withEdgeScores=True).preprocess()
						matcher = matching.PathGrowingMatcher(self.Gc[i-1], self._weightsFromADScores(ad.getEdgeScores(), self.Gc[i-1]))
					else:
						matcher = matching.PathGrowingMatcher(self.Gc[i-1])
					matcher.run()
					self.matching.append(matcher.getMatching())

					# perform coarsening
					logging.info("\t coarsening")
					coarseningAlgo = coarsening.MatchingCoarsening(self.Gc[i-1], self.matching[i-1], noSelfLoops=True)
					coarseningAlgo.run()
					self.Gc.append(coarseningAlgo.getCoarseGraph())
				else:
					raise Error("unknown aggregation scheme")

				# set node mappings
				self.down.append(coarseningAlgo.getFineToCoarseNodeMapping())
				self.up.append(coarseningAlgo.getCoarseToFineNodeMapping())

				# set node weights
				logging.info("updating node weights")
				nw = [0 for v in range(self.Gc[i].upperNodeIdBound())]
				def updateNodeWeights(v):
					nw[v] += sum(self.nodeWeights[i-1][v_] for v_ in self.up[i][v])
				self.Gc[i].forNodes(lambda v: updateNodeWeights(v))
				self.nodeWeights.append(nw)

				# stop if there's no change from level to level
				if (self.Gc[i].size() == self.Gc[i-1].size()):
					print("stopping at level ", i)
					break

	def _uncoarsen(self, i):
		""" """
		self.Gf[i] = Graph()
		Gc = self.Gf[i]
		Gf = self.Gf[i]

		for v in Gc.nodes():
			# every coarse node corresponds to a subgraph in the next finer level
			S = self.Gc[i-1].subgraphFromNodes(self.up[i][v])
			Gf.append(S)



	def _buildFineSequence(self):
		# preallocate
		self.Gf = [None for i in range(self.maxLevel + 1)]
		for i in range(self.maxLevel, -1, -1):	# count down levels
			print("level ", i)
			if i == self.maxLevel:
				self.Gf[i] = self.Gc[i]
			else:
				self._uncoarsen(i)



	def generate(self):
		self._buildCoarseSequence()
		self._buildFineSequence()
		return self.Gf[0]

	@classmethod
	def fit(cls, G):
		return cls(G)


class BTERReplicator:

	matlabname = 'matlab'
	matlabScript = """
	addpath('{0}')
	filename = 'bter_input.mat'
	load(filename)
	addpath('{1}')
	[ccd,gcc] = ccperdeg(G);
	nd = accumarray(nonzeros(sum(G,2)),1);
	nd = nd * {2};
	[E1,E2] = bter(nd,ccd);
	G_bter = bter_edges2graph(E1,E2);
	save('bter_output', 'G_bter')
	exit
	"""
	feastpackPath = "."
	workingDir = "."


	@classmethod
	def setPaths(cls, feastpackPath, workingDir="/tmp"):
		cls.feastpackPath = feastpackPath
		cls.workingDir = workingDir

	def __init__(self, G, scale=1):
		self.G = G
		self.scriptPath = os.path.join(self.workingDir, "bter_wrapper.m")
		# write MATLAB script
		with open(self.scriptPath, 'w') as matlabScriptFile:
			matlabScriptFile.write(self.matlabScript.format(self.workingDir, self.feastpackPath, scale))
		self.tempFileOut = os.path.join(self.workingDir, 'bter_output')
		self.tempFileIn = os.path.join(self.workingDir, 'bter_input.mat')

	def generate(self):
		graphio.writeMat(self.G, self.tempFileIn)
		subprocess.call([self.matlabname, '-nosplash', '-nodisplay', '-r "run(\''+self.scriptPath+'\');"'])
		G_bter = graphio.readMat(self.tempFileOut, key='G_bter')
		subprocess.call(['rm', self.tempFileOut])
		subprocess.call(['rm', self.tempFileIn])
		return G_bter

	@classmethod
	def fit(cls, G, scale=1):
		return cls(G, scale)

class MUSKETEERAdapter:

	def __init__(self, O, **params):
		self.O = nxadapter.nk2nx(O)
		defaultParams = {'node_edit_rate': [0, 0.05, 0.05, 0.05], 'edge_edit_rate': [0, 0.05, 0.05, 0.05]}
		self.params = defaultParams.copy()
		self.params.update(params)


	def generate(self):
		import sys
		currentDir = os.getcwd()
		os.chdir(self.musketeerPath)
		if sys.version_info >= (3,5):
			import importlib.util
			spec = importlib.util.spec_from_file_location("algorithms", os.path.join(self.musketeerPath, "algorithms.py"))
			musketeerModule = importlib.util.module_from_spec(spec)
			spec.loader.exec_module(musketeerModule)
		else:
			from importlib.machinery import SourceFileLoader
			musketeerModule = SourceFileLoader("algorithms", os.path.join(self.musketeerPath, "algorithms.py")).load_module()
		os.chdir(currentDir)
		R = musketeerModule.generate_graph(self.O, self.params)
		return nxadapter.nx2nk(R)

	@classmethod
	def setPaths(cls, musketeerPath):
		cls.musketeerPath = musketeerPath

	@classmethod
	def fit(cls, G, scale=1):
		params = {}
		if scale > 1:
			params = {"node_growth_rate": [0, 0.0, 0.0, scale]}
		return cls(G, **params)
