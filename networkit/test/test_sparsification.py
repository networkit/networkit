#!/usr/bin/env python3
import unittest
import os

from networkit import *

class TestSparsification(unittest.TestCase):

	def setUp(self):
		self.G = readGraph("input/jazz.graph", Format.METIS)
		self.G.indexEdges()
		self.sparsifiers = [
			sparsification.SimmelianSparsifierParametric(10),
			sparsification.SimmelianSparsifierNonParametric(),
			sparsification.QuadrilateralSimmelianSparsifier(),
			sparsification.DegreeMultiscaleSparsifier(lambda d1, d2: max(d1,d2)),
			sparsification.SimmelianMultiscaleSparsifier(),
			sparsification.LocalSimilaritySparsifier(),
			sparsification.MultiscaleSparsifier(),
			sparsification.RandomEdgeSparsifier(),
			sparsification.RandomNodeEdgeSparsifier(),
			sparsification.ForestFireSparsifier(0.6, 5.0),
			sparsification.LocalDegreeSparsifier(),
			sparsification.SCANSparsifier()
		]

	def testGetSparsifiedGraphOfSize(self):
		"""
		Checks whether the sizes of the sparsified graphs are approximately
		of the expected size. This test is supposed to verify that all
		sparsification methods can be run without errors.
		"""

		minExpectedRatio = 0.15
		targetRatio = 0.2
		maxExpectedRatio = 0.25

		for sparsifier in self.sparsifiers:
			S = sparsifier.getSparsifiedGraphOfSize(self.G, targetRatio)
			ratio = S.numberOfEdges() / self.G.numberOfEdges()
			self.assertGreaterEqual(ratio, minExpectedRatio)
			self.assertLessEqual(ratio, maxExpectedRatio)

	def testSimRankScoreRejectsInvalidSimilarityPropagationFactor(self):
		G = Graph(3, False, False)

		with self.assertRaises(ValueError):
			sparsification.SimRankScore(G, 1.1, 100, 1e-4)

		with self.assertRaises(ValueError):
			sparsification.SimRankScore(G, -0.1, 100, 1e-4)

	def testSimRankScoreRejectsZeroMaxIterations(self):
		G = Graph(3, False, False)

		with self.assertRaises(ValueError):
			sparsification.SimRankScore(G, 0.9, 0, 1e-4)

	def testSimRankScoreRejectsNegativeTolerance(self):
		G = Graph(3, False, False)

		with self.assertRaises(ValueError):
			sparsification.SimRankScore(G, 0.9, 100, -1e-4)

	def testSimRankScoreRunRequiresIndexedEdges(self):
		G = Graph(2, False, False)
		G.addEdge(0, 1)

		simrank = sparsification.SimRankScore(G)

		with self.assertRaises(RuntimeError):
			simrank.run()

	def testSimRankScoreScoresUnavailableBeforeRun(self):
		G = Graph(2, False, False)
		G.addEdge(0, 1)
		G.indexEdges()

		simrank = sparsification.SimRankScore(G)

		with self.assertRaises(RuntimeError):
			simrank.scores()

		with self.assertRaises(RuntimeError):
			simrank.score(0)

		with self.assertRaises(RuntimeError):
			simrank.score(0, 1)

	def testSimRankScoreRunAcceptsEmptyGraph(self):
		G = Graph(0, False, False)
		G.indexEdges()

		simrank = sparsification.SimRankScore(G)
		simrank.run()

		self.assertEqual(simrank.scores(), [])

	def testSimRankScoreProducesExpectedTriangleScores(self):
		G = Graph(3, False, False)
		G.addEdge(0, 1)
		G.addEdge(1, 2)
		G.addEdge(0, 2)
		G.indexEdges()

		simrank = sparsification.SimRankScore(G, 0.5, 100, 1e-12)
		simrank.run()

		scores = simrank.scores()

		self.assertEqual(len(scores), G.upperEdgeIdBound())

		for eid in range(G.upperEdgeIdBound()):
			self.assertAlmostEqual(simrank.score(eid), scores[eid])

		expected = 0.2

		for u, v in [(0, 1), (1, 2), (0, 2)]:
			eid = G.edgeId(u, v)
			self.assertAlmostEqual(simrank.score(u, v), expected, delta=1e-10)
			self.assertAlmostEqual(simrank.score(eid), expected, delta=1e-10)

	def testSimRankScoreDirectedGraphUsesInNeighbors(self):
		G = Graph(3, False, True)
		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(1, 2)
		G.indexEdges()

		simrank = sparsification.SimRankScore(G, 0.8, 1, 0.0)
		simrank.run()

		self.assertAlmostEqual(simrank.score(G.edgeId(1, 2)), 0.4, delta=1e-12)

if __name__ == "__main__":
	unittest.main()
