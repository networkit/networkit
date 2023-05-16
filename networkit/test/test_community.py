#!/usr/bin/env python3
import unittest
import os

from networkit import *


class Test_CommunityDetection(unittest.TestCase):

	def setUp(self):
		self.G = readGraph("input/PGPgiantcompo.graph",Format.METIS)

	def test_PLM(self):
		comms = community.detectCommunities(self.G, community.PLM(self.G))
		
		# each node must be assigned
		for v in self.G.iterNodes():
			self.assertTrue(comms.contains(v))

	def test_CutClustering(self):
		jazz = readGraph("input/jazz.graph",Format.METIS)

		cc = community.CutClustering(jazz,0.5076142131979697)

		comms = community.detectCommunities(jazz, cc)
		self.assertEqual(193, comms.numberOfSubsets())

		hierarchy = community.CutClustering.getClusterHierarchy(jazz)
		self.assertEqual(3, len(hierarchy))


	def test_LFM(self):
		local = scd.LFMLocal(self.G)
		lfm = community.LFM(self.G, local)
		lfm.run()
		cover = lfm.getCover()

		for u in self.G.iterNodes():
			self.assertGreaterEqual(len(cover.subsetsOf(u)), 1)

		
		
	def test_coverHubDominance(self):

		G = Graph(20)
		G.addEdge(1,2)
		G.addEdge(2,4)
		cover = Cover(20) 
		cover.setUpperBound(3)
		cover.addToSubset(0,0)
		cover.addToSubset(1,1)
		cover.addToSubset(2,2)
		
		sim=community.CoverHubDominance(G, cover)
		#sim2 = community.CoverF1Similarity(G, cover, cover)

		sim.run()
		self.assertEqual(sim.getWeightedAverage(), 1.0)
		self.assertEqual(sim.getUnweightedAverage(), 1.0)
		self.assertEqual(sim.getValue(0), 1.0)
		self.assertListEqual(sim.getValues(),[1.0, 1.0, 1.0])
		self.assertEqual(sim.getMaximumValue(), 1.0)
		self.assertEqual(sim.getMinimumValue(), 1.0)

	def test_intrapartitionDensity(self):

		G = Graph(6)
		G.addEdge(0,2)
		G.addEdge(1,2)
		G.addEdge(3,4)
	
		P = Partition(3)
		P.allToSingletons()
		
		dense=community.IntrapartitionDensity(G,P)
		dense.run()
		self.assertAlmostEqual(dense.getWeightedAverage(), 0.17, 2)
		self.assertAlmostEqual(dense.getUnweightedAverage(), 0.33, 2)
		self.assertEqual(dense.getValue(0), 1.0)
		self.assertEqual(dense.getGlobal(), 0.0)

	


if __name__ == "__main__":
	unittest.main()
