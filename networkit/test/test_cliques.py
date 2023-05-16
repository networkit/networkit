#!/usr/bin/env python3
import unittest
import os
import networkit as nk
import numpy as np
import scipy



class Test_Cliques(unittest.TestCase):

	def test_maxClique(self):
		G = nk.readGraph("input/jazz2_directed.gml",nk.Format.GML)
		G.removeSelfLoops()
		A = nk.clique.MaximalCliques(G).run().getCliques()
		
		self.assertListEqual(A,[[2], [1], [4], [3]])



if __name__ == "__main__":
	unittest.main()
