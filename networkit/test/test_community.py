import unittest
import os

from networkit import *


class Test_CommunityDetection(unittest.TestCase):

	def setUp(self):
		self.G = readGraph("input/PGPgiantcompo.graph",Format.METIS)

	def test_PLM(self):
		comms = community.detectCommunities(self.G, community.PLM(self.G))

		# each node must be assigned
		for v in self.G.nodes():
			self.assertTrue(comms.contains(v))

	def test_CutClustering(self):
		jazz = readGraph("input/jazz.graph",Format.METIS)

		cc = community.CutClustering(jazz,0.5076142131979697)

		comms = community.detectCommunities(jazz, cc)
		self.assertEqual(193, comms.numberOfSubsets())

		hierarchy = community.CutClustering.getClusterHierarchy(jazz)
		self.assertEqual(3, len(hierarchy))

