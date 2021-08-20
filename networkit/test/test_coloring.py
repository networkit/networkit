#!/usr/bin/env python3
import unittest
import os
import networkit as nk

class Test_Coloring(unittest.TestCase):

	def test_SpectralColoring(self):
		G = nk.readGraph("input/karate.graph", nk.Format.METIS)
		spCol = nk.coloring.SpectralColoring(G)

		spCol.run()

		self.assertLessEqual(len(spCol.getColoring()), G.upperNodeIdBound())

if __name__ == "__main__":
	unittest.main()
