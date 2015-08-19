import unittest
import os

from networkit import *

class Test_Sparsification(unittest.TestCase):

	def setUp(self):
		self.G = readGraph("input/jazz.graph", Format.METIS)
		self.G.indexEdges()
		self.sparsifiers = [
			sparsification.SimmelianBackboneParametric(10),
			sparsification.SimmelianBackboneNonParametric(),
			sparsification.QuadrilateralSimmelianBackbone(),
			sparsification.DegreeMultiscaleBackbone(lambda d1, d2: max(d1,d2)),
			sparsification.SimmelianMultiscaleBackbone(),
			sparsification.LocalSimilarityBackbone(),
			sparsification.MultiscaleBackbone(),
			sparsification.RandomEdgeBackbone(),
			sparsification.RandomNodeEdgeBackbone(),
			sparsification.ForestFireBackbone(0.6, 5.0),
			sparsification.LocalDegreeSparsifier(),
			sparsification.SCANBackbone()
		]

	def test_getSparsifiedGraphOfSize(self):
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
