#!/usr/bin/env python3
import numpy as np
import os
import random
import unittest

import networkit as nk

class TestComponent(unittest.TestCase):

	def setUp(self):
		self.L = nk.readGraph("input/looptest1.gml", nk.Format.GML) #without self-loops
		self.LL = nk.readGraph("input/looptest2.gml", nk.Format.GML) #with self-loops sprinkled in

	def testBiconnectedComponents(self):
		bcc = nk.components.BiconnectedComponents(self.LL)
		bcc.run()

		for component in bcc.getComponents():
			G1 = nk.graphtools.subgraphFromNodes(self.LL, component)
			def testNode(v):
				G2 = nk.Graph(G1)
				G2.removeNode(v)
				cc = nk.components.ConnectedComponents(G2)
				cc.run()
				self.assertEqual(cc.numberOfComponents(), 1)
			G1.forNodes(testNode)

	def testBiconnectedComponentsSize(self):
		G = nk.Graph(5)
		G.addEdge(0,1)
		G.addEdge(1,2)
		G.addEdge(3,4)

		bcc = nk.components.BiconnectedComponents(G)
		bcc.run()
		self.assertEqual(bcc.numberOfComponents(), 3)
		allCompSizes = bcc.getComponentSizes()
		self.assertEqual(allCompSizes[0], 2)
		self.assertEqual(allCompSizes[1], 2)
		self.assertEqual(allCompSizes[2], 2)

	def testBiconnectedComponentsNode(self):
		G = nk.Graph(5)
		G.addEdge(0,1)
		G.addEdge(1,2)
		G.addEdge(3,4)
		
		bcc = nk.components.BiconnectedComponents(G)
		bcc.run()
		self.assertEqual({1}, bcc.getComponentsOfNode(0))
		self.assertEqual({0,1}, bcc.getComponentsOfNode(1))
		self.assertEqual({0}, bcc.getComponentsOfNode(2))	
		self.assertEqual({2}, bcc.getComponentsOfNode(3))	
		self.assertEqual({2}, bcc.getComponentsOfNode(4))	
	
	def testConnectedComponents(self):
		CC = nk.components.ConnectedComponents(self.LL)
		CC.run()
		CC.componentOfNode(1)
		CC.getComponentSizes()
		CC.getPartition()
		CC.numberOfComponents()

	def testStronglyConnectedComponents(self):
		g = nk.readGraph("input/MIT8.edgelist",
				nk.Format.EdgeList, separator='\t', firstNode=0, continuous=False, directed=True)
		scc = nk.components.StronglyConnectedComponents(g)
		scc.run()
		self.assertNotEqual(scc.componentOfNode(0), None)
		nComponents = scc.numberOfComponents()
		compSizes = scc.getComponentSizes()
		self.assertEqual(nComponents, len(compSizes))

		comps = scc.getComponents()
		for idx, size in compSizes.items():
			self.assertEqual(len(comps[idx]), size)

		_=scc.getPartition()
	
	def testExtractLargestConnectedComponent(self):
		G = nk.Graph(10)
		for i in range(3):
			G.addEdge(i, i+1)

		for i in range(4, 9):
			G.addEdge(i, i+1)

		G1 = nk.components.ConnectedComponents.extractLargestConnectedComponent(G, True)
		self.assertEqual(G1.numberOfNodes(), 6)
		self.assertEqual(G1.numberOfEdges(), 5)

		G2 = nk.components.ConnectedComponents.extractLargestConnectedComponent(G, False)
		for i in range(G.numberOfNodes()):
			self.assertEqual(G2.hasNode(i), (4 <= i <= 9))

	def testParallelConnectedComponents(self):
		G = nk.Graph(10)
		for i in range(3):
			G.addEdge(i, i+1)

		for i in range(4, 9):
			G.addEdge(i, i+1)
		PCC = nk.components.ParallelConnectedComponents(G)
		PCC.run()
		#other functions are already tested extensively
		self.assertEqual(PCC.numberOfComponents(), 2)	

	def testWeaklyConnectedComponents(self):
		G = nk.Graph(10, directed=True)
		for i in range(3):
			G.addEdge(i, i+1)

		for i in range(4, 9):
			G.addEdge(i, i+1)
		WCC = nk.components.WeaklyConnectedComponents(G)
		WCC.run()
		self.assertListEqual(WCC.getComponents(),[[0, 1, 2, 3], [4, 5, 6, 7, 8, 9]])

if __name__ == "__main__":
	unittest.main()
