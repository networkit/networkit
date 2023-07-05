#!/usr/bin/env python3
import unittest
import os
import networkit as nk

class TestColoring(unittest.TestCase):

	def testSpectralColoring(self):
		G = nk.readGraph("input/karate.graph", nk.Format.METIS)
		spCol = nk.coloring.SpectralColoring(G)

		spCol.run()

		self.assertLessEqual(len(spCol.getColoring()), G.upperNodeIdBound())

if __name__ == "__main__":
	unittest.main()
