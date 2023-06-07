#!/usr/bin/env python3
import numpy as np
import os
import random
import unittest

import networkit as nk

class TestDistance(unittest.TestCase):

	def setUp(self):
		self.L = nk.readGraph("input/looptest1.gml", nk.Format.GML) #without self-loops
		self.LL = nk.readGraph("input/looptest2.gml", nk.Format.GML) #with self-loops sprinkled in

	def testAsArrayAPSP(self):
		nk.engineering.setSeed(1, True)
		random.seed(1)
		for directed in [True, False]:
			for weighted in [True, False]:
				g = nk.generators.ErdosRenyiGenerator(100, 0.15, directed).generate()
				if weighted:
					g = nk.graphtools.toWeighted(g)
					g.forEdges(lambda u, v, ew, eid: g.setWeight(u, v, random.random()))
				apsp = nk.distance.APSP(g)
				apsp.run()
				listDistances = apsp.getDistances()
				arrayDistances = apsp.getDistances(asarray=True)
				self.assertIsInstance(listDistances, list)
				self.assertIsInstance(arrayDistances, np.ndarray)
				np.testing.assert_allclose(listDistances, arrayDistances)
	
	def testAStar(self):
		# Builds a mesh graph with the given number of rows and columns
		def buildMesh(rows, cols):
			G = nk.Graph(rows * cols, False, False)
			for i in range(rows):
				for j in range(cols):
					if j < cols - 1:
						G.addEdge(i * cols + j, i * cols + j + 1)
					if i < rows - 1:
						G.addEdge(i * cols + j, (i + 1) * cols + j)
			return G

		# Test the AStar algorithm on a mesh with the given number of rows and columns
		def testMesh(rows, cols):
			G = buildMesh(rows, cols)

			# Test A* on the given source-target pair
			def testPair(s, t):

				# Some distance heuristics:

				# Always returns 0, A* degenerates to Dijkstra
				def zeroDist(u):
					return 0

				# Returns the exact distance from u to the target
				def exactDist(u):
					rowU = int(u / cols)
					colU = int(u % cols)
					rowT = int(t / cols)
					colT = int(t % cols)
					return abs(rowU - rowT) + abs(colU - colT)

				# Returns the eucledian distance from u to the target
				def eucledianDist(u):
					rowT = int(t / cols)
					colT = int(t % cols)
					rowDiff = abs(int(u / cols) - rowT)
					colDiff = abs(int(u % cols) - rowT)
					return (rowDiff**2 + colDiff**2)**.5

				# Use BFS as ground truth
				bfs = nk.distance.BFS(G, s, True, False, t).run()

				# Test A* on all the heuristics
				for heu in [zeroDist, exactDist, eucledianDist]:
					heuristics = [heu(u) for u in range(G.numberOfNodes())]
					astar = nk.distance.AStar(G, heuristics, s, t, True)
					astar.run()

					# Test distance of target
					self.assertEqual(astar.getDistance(), bfs.distance(t))

					# Test path
					path = astar.getPath()
					self.assertEqual(len(path), len(bfs.getPath(t)) - 2)
					if len(path) == 0:
						continue
					for i in range(len(path) - 1):
						self.assertTrue(G.hasEdge(path[i], path[i + 1]))

			# Iterate over all possible source-target pairs
			G.forNodePairs(testPair)

		# Test some meshes
		testMesh(10, 10)
		testMesh(21, 5)
		testMesh(9, 18)
		testMesh(7, 1)
	
	def testAdamicAdarDistance(self):
		g=nk.Graph(5, directed=False, weighted=False)
		g.addEdge(0,1)
		g.addEdge(2,3)
		g.addEdge(3,4)			
		g.indexEdges()
		adam = nk.distance.AdamicAdarDistance(g)
		adam.preprocess()
		self.assertIsInstance(adam.getAttribute(),list)
		e = nk.graphtools.randomEdge(g)
		self.assertIsInstance(adam.distance(e[0],e[1]),float)

	def testAlgebraicDistance(self):
		g=nk.Graph(5, directed=False, weighted=False)
		g.addEdge(0,1)
		g.addEdge(2,3)
		g.addEdge(3,4)
		g.indexEdges()
		alge = nk.distance.AlgebraicDistance(g, withEdgeScores=True)
		alge.preprocess()
		self.assertIsInstance(alge.getEdgeScores(),list)
		e = nk.graphtools.randomEdge(g)
		self.assertIsInstance(alge.distance(e[0],e[1]),float)			
		
	def testCommuteTimeDistance(self):
		g=nk.Graph(5, directed=False, weighted=False)
		g.addEdge(0,1)
		g.addEdge(2,3)
		g.addEdge(3,4)
		g.indexEdges()
		ctd = nk.distance.CommuteTimeDistance(g)
		ctd.runApproximation()
		ctd.runParallelApproximation()
		ctd.runSinglePair(nk.graphtools.randomNode(g),nk.graphtools.randomNode(g))
		ctd.runSingleSource(nk.graphtools.randomNode(g))
		e = nk.graphtools.randomEdge(g)
		self.assertIsInstance(ctd.distance(e[0],e[1]),float)	
	
	def testDiameterRange(self):
		D = nk.distance.Diameter(self.LL, nk.distance.DiameterAlgo.ESTIMATED_RANGE, error = 0.1)
		D.run()
		self.assertIsInstance(D.getDiameter(), tuple)

	def testDiameterSamples(self):		
		D = nk.distance.Diameter(self.LL, nk.distance.DiameterAlgo.ESTIMATED_SAMPLES, nSamples = 5)
		D.run()
		self.assertIsInstance(D.getDiameter(), tuple)

	def testDiameterRangeExact(self):		
		D = nk.distance.Diameter(self.LL, nk.distance.DiameterAlgo.EXACT)
		D.run()
		self.assertTupleEqual(D.getDiameter(), (5,0))

	def testEccentricity(self):
		E = nk.distance.Eccentricity()
		E.getValue(self.LL, 0)

	def testEffectiveDiameter(self):
		algo = nk.distance.EffectiveDiameter(self.L)
		algo.run()
		self.assertEqual(algo.getEffectiveDiameter(), 4.0)
		algo = nk.distance.EffectiveDiameter(self.LL)
		algo.run()
		self.assertEqual(algo.getEffectiveDiameter(), 4.0)

	def testApproxEffectiveDiameter(self):
		algo = nk.distance.EffectiveDiameterApproximation(self.L)
		algo.run()
		self.assertAlmostEqual(algo.getEffectiveDiameter(), 3.55, delta=1)
		algo = nk.distance.EffectiveDiameterApproximation(self.LL)
		algo.run()
		self.assertAlmostEqual(algo.getEffectiveDiameter(), 3.55, delta=1)

	def testApproxHopPlot(self):
		algo = nk.distance.HopPlotApproximation(self.L)
		algo.run()
		self.assertAlmostEqual(len(algo.getHopPlot()), 6, delta=2)
		algo = nk.distance.HopPlotApproximation(self.LL)
		algo.run()
		self.assertAlmostEqual(len(algo.getHopPlot()), 6, delta=2)

	def testNeighborhoodFunction(self):
		algo = nk.distance.NeighborhoodFunction(self.L)
		algo.run()
		self.assertListEqual(algo.getNeighborhoodFunction(), [24, 44, 60, 70, 72])
		algo = nk.distance.NeighborhoodFunction(self.LL)
		algo.run()
		self.assertListEqual(algo.getNeighborhoodFunction(), [24, 44, 60, 70, 72])

	def testApproxNeighborhoodFunction(self):
		algo = nk.distance.NeighborhoodFunctionApproximation(self.L)
		algo.run()
		self.assertEqual(len(algo.getNeighborhoodFunction()), 5)
		algo = nk.distance.NeighborhoodFunctionApproximation(self.LL)
		algo.run()
		self.assertEqual(len(algo.getNeighborhoodFunction()), 5)

	def genERGraphs(self, n = 100, p = 0.15, seed = 42):
		nk.engineering.setSeed(seed, True)
		random.seed(seed)
		for directed in [True, False]:
			for weighted in [True, False]:
				g = nk.generators.ErdosRenyiGenerator(100, 0.15, directed).generate()
				if weighted:
					g = nk.graphtools.toWeighted(g)
					g.forEdges(lambda u, v, ew, eid: g.setWeight(u, v, random.random()))
				yield g

	def testJaccardSimilarity(self):	
		for g in self.genERGraphs():
			if g.isDirected() and not g.isWeighted():
				g.indexEdges()
				jaccSim = nk.distance.JaccardSimilarityAttributizer(g,[1])
				#no complete list comparison, length > 1000
				res = jaccSim.getAttribute()[0:2]
				self.assertAlmostEqual(res[0], 0.038, delta=0.01)
				self.assertEqual(res[1], 0.0)
	
	def testSPSP(self):
		for g in self.genERGraphs():
			for nSources in [1, 10, 50]:
				sources = nk.graphtools.randomNodes(g, nSources)
				spsp = nk.distance.SPSP(g, sources)
				spsp.run()
				dists = spsp.getDistances()
				self.assertEqual(len(dists), nSources)
				for distList in dists:
					self.assertEqual(len(distList), g.numberOfNodes())

				for nTargets in [1, 50, 100]:
					targets = nk.graphtools.randomNodes(g, nTargets)
					spsp = nk.distance.SPSP(g, sources, targets)
					spsp.run()
					dists = spsp.getDistances()
					self.assertEqual(len(dists), nSources)
					for distList in dists:
						self.assertEqual(len(distList), nTargets)

	def testMultiTargetSTSP(self):
		for g in self.genERGraphs():
			source = nk.graphtools.randomNode(g)
			for nTargets in [1, 10, 50]:
				targets = nk.graphtools.randomNodes(g, nTargets)
				algo = None
				if g.isWeighted():
					algo = nk.distance.MultiTargetDijkstra(g, source, targets)
				else:
					algo = nk.distance.MultiTargetBFS(g, source, targets)
				algo.setSource(source)
				algo.setTarget(targets[nTargets-1])
				algo.run()
				self.assertLessEqual(len(algo.getPredecessors()),g.numberOfNodes())

	def testSSSP(self):
		for g in self.genERGraphs():
			source = nk.graphtools.randomNode(g)
			for nTargets in [1, 10, 50]:
				targets = nk.graphtools.randomNodes(g, nTargets)
				algo = None
				if g.isWeighted():
					algo = nk.distance.Dijkstra(g, source, targets, storeNodesSortedByDistance = True)
				else:
					algo = nk.distance.BFS(g, source, targets)
				algo.setSource(source)	
				algo.setTarget(targets[nTargets-1])
				algo.run()
				self.assertLessEqual(len(algo.getPredecessors(source)),g.numberOfNodes())
				self.assertLessEqual(len(targets),len(algo.getDistances()))
				self.assertLessEqual(len(algo.getPaths(targets[nTargets-1])),g.numberOfNodes())
				if g.isWeighted():
					self.assertLessEqual(len(algo.getNodesSortedByDistance()), g.numberOfNodes())
				self.assertLessEqual(1.0, algo.numberOfPaths(targets[nTargets-1]))

	def testDynSSSP(self):
		for g in self.genERGraphs():
			source = nk.graphtools.randomNode(g)
			for nTargets in [1, 10, 50]:
				targets = nk.graphtools.randomNodes(g, nTargets)
				algoDyn = None
				if g.isWeighted():
					algoDyn = nk.distance.DynDijkstra(g, source)
				else:
					algoDyn = nk.distance.DynBFS(g, source)
				algoDyn.run()
				if not g.isWeighted():
					batch1 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 3, 7, 1.0)
					batch2 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 1, 5, 1.0)
					batch = [batch1,batch2]
					algoDyn.updateBatch(batch)
					algoDyn.run()
					self.assertIsInstance(algoDyn.modified(), bool)

	def testPrunedLandmarkLabeling(self):
		for g in self.genERGraphs():
			pll = nk.distance.PrunedLandmarkLabeling(g)
			pll.run()

			if g.isWeighted():
				g = nk.graphtools.toUnweighted(g)
			apsp = nk.distance.APSP(g)
			apsp.run()

			for u in g.iterNodes():
				for v in g.iterNodes():
					if apsp.getDistance(u, v) > g.numberOfNodes():
						self.assertEqual(pll.query(u, v), nk.none)
					else:
						self.assertEqual(pll.query(u, v), int(apsp.getDistance(u, v)))

	def testVolume(self):
		vol = nk.distance.Volume()
		self.assertEqual(vol.volume(self.L, 0.5), 1.0)

	def testWithNonConsecutiveNodeLabelsAPSP(self):
		nk.engineering.setSeed(1, True)
		random.seed(1)
		g = nk.generators.ErdosRenyiGenerator(100, 0.15).generate()
		g.removeNode(0)
		nk.distance.APSP(g).run()		

if __name__ == "__main__":
	unittest.main()
