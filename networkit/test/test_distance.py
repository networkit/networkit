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
	
	def testAlgebraicDistance(self):
		nk.setSeed(42, True)
		g=nk.Graph(5, directed=False, weighted=False)
		g.addEdge(0,1)
		g.addEdge(2,3)
		g.addEdge(3,4)
		g.indexEdges()
		alge = nk.distance.AlgebraicDistance(g, withEdgeScores=True)
		alge.preprocess()
		res = [0.0, 2.053536340973494e-09, 2.0535364519957966e-09]
		for i in range(len(res)):
			self.assertAlmostEqual(alge.getEdgeScores()[i], res[i], delta=0.01)
		self.assertAlmostEqual(alge.distance(2,4), 4.107072792969291e-09, delta=0.01)			
		
	def testCommuteTimeDistanceApprox(self):
		g=nk.Graph(5, directed=False, weighted=False)
		g.addEdge(0,1)
		g.addEdge(4,1)
		g.addEdge(2,3)
		g.addEdge(3,4)
		g.indexEdges()
		ctd = nk.distance.CommuteTimeDistance(g)
		ctd.runApproximation()
		self.assertAlmostEqual(ctd.distance(2,4), 4.0990815830507925, delta=0.5)

	def testCommuteTimeDistanceParallelApprox(self):
		g=nk.Graph(5, directed=False, weighted=False)
		g.addEdge(0,1)
		g.addEdge(4,1)
		g.addEdge(2,3)
		g.addEdge(3,4)
		g.indexEdges()
		ctd = nk.distance.CommuteTimeDistance(g)
		ctd.runParallelApproximation()
		self.assertAlmostEqual(ctd.distance(1,4), 2.8284089376767025, delta=0.5)

	def testCommuteTimeDistanceSinglePair(self):
		g=nk.Graph(5, directed=False, weighted=False)
		g.addEdge(0,1)
		g.addEdge(4,1)
		g.addEdge(2,3)
		g.addEdge(3,4)
		g.indexEdges()
		ctd = nk.distance.CommuteTimeDistance(g).run()
		ctd.runSinglePair(2,4)
		self.assertAlmostEqual(ctd.distance(2,4), 4.0, delta=0.5)

	def testCommuteTimeDistanceSingleSource(self):
		g=nk.Graph(5, directed=False, weighted=False)
		g.addEdge(0,1)
		g.addEdge(4,1)
		g.addEdge(2,3)
		g.addEdge(3,4)
		g.indexEdges()
		ctd = nk.distance.CommuteTimeDistance(g).run()
		ctd.runSingleSource(3)
		self.assertAlmostEqual(ctd.distance(1,3), 4.0, delta=0.5)

	def testDiameterRange(self):
		D = nk.distance.Diameter(self.LL, nk.distance.DiameterAlgo.ESTIMATED_RANGE, error = 0.01)
		D.run()
		self.assertTupleEqual(D.getDiameter(), (5,5))

	def testDiameterSamples(self):		
		D = nk.distance.Diameter(self.LL, nk.distance.DiameterAlgo.ESTIMATED_SAMPLES, nSamples = 100)
		D.run()
		self.assertTupleEqual(D.getDiameter(), (10,0))


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

	def testApproxEffectiveDiameter(self):
		algo = nk.distance.EffectiveDiameterApproximation(self.L)
		algo.run()
		self.assertAlmostEqual(algo.getEffectiveDiameter(), 3.55, delta=1)

	def testApproxHopPlot(self):
		algo = nk.distance.HopPlotApproximation(self.L)
		algo.run()
		self.assertAlmostEqual(len(algo.getHopPlot()), 6, delta=2)

	def testNeighborhoodFunction(self):
		algo = nk.distance.NeighborhoodFunction(self.L)
		algo.run()
		self.assertListEqual(algo.getNeighborhoodFunction(), [24, 44, 60, 70, 72])

	def testApproxNeighborhoodFunction(self):
		algo = nk.distance.NeighborhoodFunctionApproximation(self.L)
		algo.run()
		self.assertEqual(len(algo.getNeighborhoodFunction()), 5)

	def testApproxNeighborhoodFunctionHeuristic(self):
		algo = nk.distance.NeighborhoodFunctionHeuristic(self.L)
		algo.run()
		self.assertEqual(algo.getNeighborhoodFunction(), [24, 42, 62, 73, 72])
	
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

	def testAsArraySPSP(self):
		nk.engineering.setSeed(1, True)
		random.seed(1)
		for directed in [True, False]:
			for weighted in [True, False]:
				g = nk.generators.ErdosRenyiGenerator(100, 0.15, directed).generate()
				if weighted:
					g = nk.graphtools.toWeighted(g)
					g.forEdges(lambda u, v, ew, eid: g.setWeight(u, v, random.random()))
				spsp = nk.distance.SPSP(g, [0,1,2,3,4])
				spsp.run()
				listDistances = spsp.getDistances()
				arrayDistances = spsp.getDistances(asarray=True)
				self.assertIsInstance(listDistances, list)
				self.assertIsInstance(arrayDistances, np.ndarray)
				np.testing.assert_allclose(listDistances, arrayDistances)
	
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
				if g.isWeighted():
					algoDyn = nk.distance.DynDijkstra(g, source)
				else:
					algoDyn = nk.distance.DynBFS(g, source)
				algoDyn.run()
				if not g.isWeighted():
					batch = [nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 3, 7, 1.0)]
					batch.append(nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 1, 5, 1.0))
					algoDyn.updateBatch(batch)
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

	def testDynPrunedLandmarkLabeling(self):
		# 0       3
		#  \     / \
		#   1---2   4

		g = nk.Graph(5, False, False)
		for e in [[0, 1], [2, 3], [3, 4]]:
			g.addEdge(e[0], e[1])

		pll = nk.distance.DynPrunedLandmarkLabeling(g)
		pll.run()

		self.assertEqual(pll.query(0, 1), 1)
		self.assertEqual(pll.query(0, 2), nk.none)

		g.addEdge(1, 2)
		pll.update(nk.dynamics.GraphEvent(
			nk.dynamics.GraphEventType.EDGE_ADDITION, 1, 2, 1))

		self.assertEqual(pll.query(0, 1), 1)
		self.assertEqual(pll.query(2, 4), 2)

	def testAdamicAdarDistance(self):
		self.L.indexEdges()
		adam = nk.distance.AdamicAdarDistance(self.L)
		adam.preprocess()
		self.assertEqual(len(adam.getAttribute()), 12)
		self.assertAlmostEqual(adam.distance(2,3), 0.6931471805599453, delta=0.1)	
	
	def testJaccardDistance(self):	
		jaccDis = nk.distance.JaccardDistance(self.L,[1,2,3])
		self.assertIsInstance(jaccDis.getAttribute(), list)	
	
	def testJaccardSimilarity(self):	
		G=nk.Graph(5)
		G.addEdge(0,1)
		G.addEdge(2,1)
		G.addEdge(3,4)
		G.indexEdges()
		jaccSim = nk.distance.JaccardSimilarityAttributizer(G,[1,2,3])
		self.assertListEqual(jaccSim.getAttribute(), [0.5, 2.0, 0.0])	
	
	def testReverseBFS(self):
		rBFS = nk.distance.ReverseBFS(self.L, 0)
		rBFS.run()
		self.assertEqual(rBFS.numberOfPaths(8), 1.0)
	
	def testVolume(self):
		vol = nk.distance.Volume()
		#single instance
		self.assertEqual(vol.volume(self.L, 0.5), 1.0)
		#list test
		self.assertListEqual(vol.volume(self.L, [0.1, 0.2, 0.8]), [1.0, 1.0, 1.0])
		#assert correct types
		self.assertRaises(Exception, vol.volume(self.L, [0.1, 0.2, "wrongType"]))

if __name__ == "__main__":
	unittest.main()
