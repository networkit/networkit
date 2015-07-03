import unittest
import os

from networkit import *

class Test_Sparsification(unittest.TestCase):

	def setUp(self):
		self.G = readGraph("input/MIT8.edgelist", Format.EdgeListTabZero)
		self.G.indexEdges()
		self.sparsifiers = [
			sparsification.SimmelianBackboneParametric(10),
			sparsification.SimmelianBackboneNonParametric(),
			sparsification.SimmelianMultiscaleBackbone(),
			sparsification.LocalSimilarityBackbone(),
			sparsification.MultiscaleBackbone(),
			sparsification.RandomEdgeBackbone(),
			sparsification.ForestFireBackbone(0.6, 5.0),
			sparsification.LocalDegreeBackbone()
		]

	def test_getSparsifiedGraphOfSize(self):
		#Verify that the size of the sparsified graphs is approximately the expected size.

		minExpectedRatio = 0.15
		targetRatio = 0.2
		maxExpectedRatio = 0.25

		for sparsifier in self.sparsifiers:
			S = sparsifier.getSparsifiedGraphOfSize(self.G, targetRatio)
			ratio = S.numberOfEdges() / self.G.numberOfEdges()
			self.assertTrue(ratio >= minExpectedRatio)
			self.assertTrue(ratio <= maxExpectedRatio)
