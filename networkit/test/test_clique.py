#!/usr/bin/env python3
import unittest
import os
import networkit as nk
import numpy as np
import scipy

class TestCliques(unittest.TestCase):

	
	def testMaxClique(self):
		G = nk.readGraph("input/jazz2_directed.gml",nk.Format.GML)
		G.removeSelfLoops()
		#without callback, testing maxOnly, expecting single Clique
		A1 = nk.clique.MaximalCliques(G, maximumOnly=True).run().getCliques()
		self.assertListEqual(A1, [[2]])
		
	def testMaxCliqueWithCallback(self):
		G = nk.readGraph("input/karate.graph",nk.Format.METIS)
		G.removeSelfLoops()
		
		def checkCliqueNodes(maxClique):
			for u in maxClique:
				self.assertTrue(G.hasNode(u))
				
		#with callback, testing allMaxCliques
		nk.clique.MaximalCliques(G, maximumOnly=False, callback=checkCliqueNodes).run()

if __name__ == "__main__":
	unittest.main()